"""Offline exact-key scoring pipeline over validated M2 observations."""

import json
import shutil
from collections import defaultdict
from decimal import Decimal
from hashlib import sha256
from pathlib import Path
from typing import Any

import polars as pl

from .bundle import validate_polygenic_score_bundle
from .models import Evaluability, ScoreConfig, ScoreRunResult, decimal_string


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _dump(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def _validate_m2(root: Path) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    manifest = json.loads((root / "manifest.json").read_text())
    for name, digest in manifest.get("artifacts", {}).items():
        if not (root / name).is_file() or _hash(root / name) != digest:
            raise ValueError(f"M2 artifact checksum mismatch: {name}")
    variants = pl.read_parquet(root / "variants.parquet").to_dicts()
    genotypes = pl.read_parquet(root / "canonical_genotypes.parquet").to_dicts()
    if any(x["normalization_run_id"] != manifest["run_id"] for x in genotypes):
        raise ValueError("M2 genotype lineage mismatch")
    return manifest, variants, genotypes


def calculate_polygenic_scores(
    normalization_directory: Path,
    score_bundle_directory: Path,
    output_directory: Path,
    config: ScoreConfig,
    imputation_directory: Path | None = None,
) -> ScoreRunResult:
    """Calculate selected scores without network access or analysis-time harmonization."""
    bundle = validate_polygenic_score_bundle(score_bundle_directory)
    unknown = set(config.pgs_ids) - set(bundle.model_ids)
    if unknown:
        raise ValueError(f"requested PGS IDs absent from checked bundle: {sorted(unknown)}")
    if imputation_directory is not None or config.genotype_source_policy != "observed_only":
        raise ValueError(
            "M6 dosage input remains gated until the completed-run validator is available"
        )
    manifest, variants, genotypes = _validate_m2(normalization_directory)
    by_variant = {v["variant_id"]: v for v in variants}
    calls: dict[tuple[str, str, int, str, str], set[tuple[str, ...]]] = defaultdict(set)
    for gt in genotypes:
        v = by_variant.get(gt["variant_id"])
        if v and v["assembly"] == "GRCh38" and gt["ploidy"] == 2 and gt["call_status"] == "called":
            calls[
                (v["assembly"], v["chromosome"], v["position"], v["reference"], v["alternate"])
            ].add(tuple(sorted(gt["alleles"])))
    metadata = json.loads((bundle.directory / "model_metadata.json").read_text())
    model_meta = {x["pgs_id"]: x for x in metadata["models"]}
    model_rows = [
        x
        for x in pl.read_parquet(bundle.directory / "model_variants.parquet").to_dicts()
        if x["pgs_id"] in config.pgs_ids
    ]
    alignment, contributions, exclusions, results, references = [], [], [], [], []
    for pgs_id in config.pgs_ids:
        rows = [x for x in model_rows if x["pgs_id"] == pgs_id]
        total_abs, matched_abs, score = Decimal(0), Decimal(0), Decimal(0)
        contributed = 0
        for row in rows:
            weight = Decimal(row["effect_weight"])
            total_abs += abs(weight)
            key = (
                row["assembly"],
                row["chromosome"],
                row["position"],
                row["reference"],
                row["alternate"],
            )
            observed = calls.get(key, set())
            if not observed:
                reason = "TARGET_VARIANT_MISSING"
            elif len(observed) != 1:
                reason = "DUPLICATE_OBSERVATION_CONFLICT"
            else:
                alleles = next(iter(observed))
                if any(a not in (row["reference"], row["alternate"]) for a in alleles):
                    reason = "TARGET_ALLELE_INCOMPATIBLE"
                else:
                    alt = Decimal(sum(a == row["alternate"] for a in alleles))
                    dosage = alt if row["effect_allele"] == row["alternate"] else Decimal(2) - alt
                    value = weight * dosage
                    score += value
                    matched_abs += abs(weight)
                    contributed += 1
                    alignment.append(
                        {
                            "pgs_id": pgs_id,
                            "variant_key": "|".join(map(str, key)),
                            "status": "matched",
                            "reason": None,
                            "input_source": "observed",
                        }
                    )
                    contributions.append(
                        {
                            "pgs_id": pgs_id,
                            "variant_key": "|".join(map(str, key)),
                            "effect_allele": row["effect_allele"],
                            "effect_weight_original": row["effect_weight"],
                            "effect_weight": decimal_string(weight),
                            "effect_dosage": decimal_string(dosage),
                            "input_source": "observed",
                            "contribution": decimal_string(value),
                            "m2_run_id": manifest["run_id"],
                        }
                    )
                    continue
            alignment.append(
                {
                    "pgs_id": pgs_id,
                    "variant_key": "|".join(map(str, key)),
                    "status": "excluded",
                    "reason": reason,
                    "input_source": None,
                }
            )
            exclusions.append(
                {"pgs_id": pgs_id, "variant_key": "|".join(map(str, key)), "reason": reason}
            )
        complete = contributed == len(rows) and len(rows) > 0
        policy = model_meta[pgs_id]["completeness_policy"]
        status = (
            Evaluability.EVALUABLE
            if complete and policy == "all_supported_variants_required"
            else Evaluability.PARTIAL_NOT_EVALUABLE
            if contributed
            else Evaluability.NOT_EVALUABLE
        )
        result = {
            "pgs_id": pgs_id,
            "model_version": model_meta[pgs_id]["version"],
            "status": status.value,
            "raw_partial_score": decimal_string(score),
            "supported_marker_count": len(rows),
            "contributed_marker_count": contributed,
            "missing_or_excluded_marker_count": len(rows) - contributed,
            "observed_contributed_count": contributed,
            "imputed_contributed_count": 0,
            "marker_coverage": contributed / len(rows) if rows else 0.0,
            "absolute_weight_coverage": float(matched_abs / total_abs) if total_abs else 0.0,
        }
        results.append(result)
        references.append(
            {
                "pgs_id": pgs_id,
                **(
                    {"status": "not_evaluable", "reason": "MODEL_NOT_EVALUABLE"}
                    if status != Evaluability.EVALUABLE
                    else {"status": "not_evaluable", "reason": "REFERENCE_NOT_SELECTED"}
                ),
            }
        )
    identity = {
        "m2_run_id": manifest["run_id"],
        "bundle_id": bundle.bundle_id,
        "config": config.model_dump(mode="json"),
        "results": results,
    }
    run_id = "m7-" + sha256(_dump(identity)).hexdigest()
    if output_directory.exists():
        raise FileExistsError("M7 output already exists")
    temp = output_directory.with_name(f".{output_directory.name}.tmp")
    temp.mkdir(parents=True)
    try:
        table_data = {
            "score_models.parquet": [model_meta[x] for x in config.pgs_ids],
            "score_variant_alignment.parquet": alignment,
            "score_contributions.parquet": contributions,
            "score_exclusions.parquet": exclusions,
            "polygenic_score_results.parquet": results,
            "reference_comparisons.parquet": references,
        }
        artifacts = {}
        for name, rows in table_data.items():
            pl.DataFrame(rows).write_parquet(temp / name)
            artifacts[name] = _hash(temp / name)
        qc = {
            "schema": "genome-evidence-pgs-qc/v1",
            "model_count": len(results),
            "statuses": {x["pgs_id"]: x["status"] for x in results},
        }
        for name, data in (
            ("metadata.json", _dump(identity)),
            ("polygenic_score_qc.json", _dump(qc)),
            (
                "report.md",
                (
                    (
                        "# Polygenic score calculation\n\n"
                        "Aggregate-only report. A raw or partial score "
                        "is not a probability, diagnosis, absolute risk, treatment recommendation, "
                        "or portable percentile. Performance and calibration can vary "
                        "across cohorts.\n\n"
                    )
                    + "\n".join(
                        f"- {x['pgs_id']}: {x['status']}; "
                        f"marker coverage {x['marker_coverage']:.3f}"
                        for x in results
                    )
                    + "\n"
                ).encode(),
            ),
        ):
            (temp / name).write_bytes(data)
            artifacts[name] = _hash(temp / name)
        (temp / "manifest.json").write_bytes(
            _dump(
                {
                    "schema": "genome-evidence-m7-manifest/v1",
                    "run_id": run_id,
                    "artifacts": artifacts,
                }
            )
        )
        artifacts["manifest.json"] = _hash(temp / "manifest.json")
        completion = {
            "schema": "genome-evidence-m7-completion/v1",
            "run_id": run_id,
            "artifact_hashes": artifacts,
        }
        (temp / "COMPLETED.json").write_bytes(_dump(completion))
        temp.rename(output_directory)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    return ScoreRunResult(
        run_id=run_id,
        output_directory=output_directory,
        model_count=len(results),
        statuses={x["pgs_id"]: Evaluability(x["status"]) for x in results},
    )
