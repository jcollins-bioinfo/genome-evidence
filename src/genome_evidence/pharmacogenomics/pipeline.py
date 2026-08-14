"""Offline M2-to-M8 orchestration and immutable aggregate-first publication."""

import json
import shutil
import subprocess
from collections import Counter
from hashlib import sha256
from pathlib import Path
from typing import Any

import polars as pl

from genome_evidence import __version__

from .bundle import validate_pharmacogenomics_bundle
from .matching import enumerate_diplotypes, observed_locus_evidence
from .models import MatchingLimits, PharmacogenomicsResult

RUN_SCHEMA = "genome-evidence-m8-run/v1"
COMPLETION_SCHEMA = "genome-evidence-m8-completion/v1"


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _dump(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, default=str) + "\n").encode()


def _validate_m2(run: Path) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    try:
        manifest = json.loads((run / "manifest.json").read_text())
    except Exception as error:
        raise ValueError("normalization run is not a complete M2 run") from error
    required = {"variants.parquet", "canonical_genotypes.parquet", "normalization_metadata.json"}
    if manifest.get("schema_version") != 1 or not required <= set(manifest.get("artifacts", {})):
        raise ValueError("normalization run has an unknown or incomplete schema")
    for name, digest in manifest["artifacts"].items():
        path = run / name
        if not path.is_file() or path.is_symlink() or _hash(path) != digest:
            raise ValueError("normalization run artifact integrity failure")
    metadata = json.loads((run / "normalization_metadata.json").read_text())
    if (
        metadata.get("run_id") != manifest.get("run_id")
        or metadata.get("target_assembly") != "GRCh38"
    ):
        raise ValueError("normalization run identity or assembly is incompatible")
    variants = pl.read_parquet(run / "variants.parquet").to_dicts()
    genotypes = pl.read_parquet(run / "canonical_genotypes.parquet").to_dicts()
    if any(row["normalization_run_id"] != manifest["run_id"] for row in genotypes):
        raise ValueError("canonical genotype lineage does not match M2 run")
    return manifest, variants, genotypes


def infer_pharmacogenomics(
    normalization_run: Path,
    bundle: Path,
    output: Path,
    genes: tuple[str, ...],
    limits: MatchingLimits | None = None,
) -> PharmacogenomicsResult:
    """Infer conservative candidate evidence from validated M2 observations.

    Analysis is local and performs no acquisition or network access. Selected genes
    are mandatory; missing positions are not filled. Publication stages files in a
    sibling directory, refuses overwrite, and writes ``COMPLETED.json`` last.
    """
    if not genes or len(genes) != len(set(genes)):
        raise ValueError("select one or more unique pharmacogenes explicitly")
    validation = validate_pharmacogenomics_bundle(bundle)
    unknown = sorted(set(genes) - set(validation.genes))
    if unknown:
        raise ValueError("selected gene is not declared by the validated bundle")
    m2_manifest, variants, genotypes = _validate_m2(normalization_run)
    chosen_limits = limits or MatchingLimits()
    results = tuple(
        enumerate_diplotypes(
            validation.directory,
            gene,
            observed_locus_evidence(validation.directory, variants, genotypes, gene),
            chosen_limits,
        )
        for gene in sorted(genes)
    )
    scientific = {
        "m2_run_id": m2_manifest["run_id"],
        "m2_manifest_sha256": _hash(normalization_run / "manifest.json"),
        "bundle_id": validation.bundle_id,
        "bundle_hash": validation.bundle_hash,
        "genes": sorted(genes),
        "limits": chosen_limits.model_dump(),
        "algorithm": "m8-exact-unphased-v1",
    }
    run_id = "m8-" + sha256(_dump(scientific)).hexdigest()
    if output.exists():
        raise FileExistsError("pharmacogenomics output already exists")
    stage = output.with_name(f".{output.name}.{run_id}.tmp")
    stage.mkdir(parents=True)
    try:
        tables: dict[str, list[dict[str, Any]]] = {
            "gene_evaluability.parquet": [
                {
                    "gene_id": item.gene_id,
                    "outcome": item.outcome.value,
                    "candidate_count": len(item.candidates),
                }
                for item in results
            ],
            "locus_evidence.parquet": [
                {"gene_id": item.gene_id, **row.model_dump(mode="json")}
                for item in results
                for row in item.locus_evidence
            ],
            "allele_candidates.parquet": [],
            "diplotype_candidates.parquet": [
                row.model_dump(mode="json") for item in results for row in item.candidates
            ],
            "candidate_exclusions.parquet": [],
            "allele_function_evidence.parquet": [],
            "phenotype_evidence.parquet": [],
            "guideline_evidence_links.parquet": [],
        }
        schemas = {
            "allele_candidates.parquet": {"gene_id": pl.String, "allele_id": pl.String},
            "candidate_exclusions.parquet": {"gene_id": pl.String, "reason": pl.String},
            "allele_function_evidence.parquet": {
                "candidate_id": pl.String,
                "source_term": pl.String,
            },
            "phenotype_evidence.parquet": {
                "gene_id": pl.String,
                "status": pl.String,
                "source_term": pl.String,
            },
            "guideline_evidence_links.parquet": {
                "gene_id": pl.String,
                "guideline_id": pl.String,
                "url": pl.String,
            },
        }
        hashes: dict[str, dict[str, object]] = {}
        for name, rows in tables.items():
            frame = pl.DataFrame(rows) if rows else pl.DataFrame(schema=schemas.get(name, {}))
            frame.write_parquet(stage / name)
            hashes[name] = {
                "sha256": _hash(stage / name),
                "byte_size": (stage / name).stat().st_size,
                "row_count": frame.height,
            }
        counts = dict(sorted(Counter(item.outcome.value for item in results).items()))
        metadata = {
            "schema": RUN_SCHEMA,
            "run_id": run_id,
            "package_version": __version__,
            "git_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], capture_output=True, text=True
            ).stdout.strip(),
            **scientific,
        }
        qc = {
            "selected_gene_count": len(genes),
            "gene_outcomes": counts,
            "candidate_count": sum(len(item.candidates) for item in results),
        }
        report = (
            "# Research-only pharmacogenomic candidate evidence\n\n"
            "> Not a clinical test, prescribing engine, medication recommendation, "
            "or replacement for a validated laboratory assay.\n\n"
            "## Aggregate evaluability\n"
            + "".join(f"- {key}: {value}\n" for key, value in counts.items())
            + "\n## Limitations\nMissing loci never become reference or `*1`. "
            "Unphased ambiguity is retained. Structural alleles, CYP2D6, HLA, "
            "mitochondrial and other special methods are unsupported. Medication "
            "decisions require an appropriate validated assay and qualified clinician/"
            "pharmacist review of the current guideline and full clinical context.\n"
        )
        for name, content in (
            ("pharmacogenomics_metadata.json", _dump(metadata)),
            ("pharmacogenomics_qc.json", _dump(qc)),
            ("pharmacogenomics_report.md", report.encode()),
        ):
            (stage / name).write_bytes(content)
            hashes[name] = {
                "sha256": _hash(stage / name),
                "byte_size": (stage / name).stat().st_size,
            }
        manifest = {"schema": RUN_SCHEMA, "run_id": run_id, "artifacts": hashes}
        (stage / "manifest.json").write_bytes(_dump(manifest))
        completion = {
            "schema": COMPLETION_SCHEMA,
            "run_id": run_id,
            "manifest_sha256": _hash(stage / "manifest.json"),
        }
        (stage / "COMPLETED.json").write_bytes(_dump(completion))
        stage.rename(output)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return PharmacogenomicsResult(run_id=run_id, output_directory=output, gene_results=results)
