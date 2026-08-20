"""M9 orchestration and atomic, checksum-declared artifact publication."""

import json
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from typing import Any

import polars as pl

from genome_evidence import __version__

from .inputs import load_m2_input
from .models import (
    ALGORITHM_VERSION,
    RUN_SCHEMA,
    CompatibilityStatus,
    FamilyAnalysisResult,
)
from .pedigree import load_pedigree
from .segregation import evaluate_site


def _dump(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, default=str) + "\n").encode()


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _table(path: Path, rows: list[dict[str, Any]], schema: dict[str, type[pl.DataType]]) -> None:
    frame = pl.DataFrame(rows) if rows else pl.DataFrame(schema=schema)
    frame.write_parquet(path)


def analyze_family(pedigree: Path, output: Path) -> FamilyAnalysisResult:
    """Analyze declared trios/duos using exact observed M2 canonical genotype identity.

    Publication contains pseudonymous lineage and site evidence and is therefore
    private. Default CLI output reports aggregates only. ``COMPLETED.json`` is
    written last after every declared artifact has been checksummed.
    """
    descriptor = load_pedigree(pedigree)
    inputs = tuple(load_m2_input(member, pedigree.parent) for member in descriptor.members)
    assemblies = {item.metadata.get("target_assembly") for item in inputs}
    algorithms = {item.metadata.get("algorithm") for item in inputs}
    if assemblies != {"GRCh38"} or algorithms != {"m2-snv-1"}:
        raise ValueError("family inputs use mixed assembly or canonical representation")
    scientific = {
        "pedigree_sha256": sha256(
            _dump(descriptor.model_dump(mode="json", by_alias=True))
        ).hexdigest(),
        "pedigree_schema": descriptor.schema_id,
        "family_id": descriptor.family_id,
        "inputs": sorted(
            (item.member.member_id, item.run_id, item.manifest_sha256) for item in inputs
        ),
        "relationships": [edge.model_dump(mode="json") for edge in descriptor.relationships],
        "algorithm": ALGORITHM_VERSION,
    }
    run_id = "m9-" + sha256(_dump(scientific)).hexdigest()
    if output.exists():
        raise FileExistsError("family-analysis output already exists")
    by_member = {item.member.member_id: item for item in inputs}
    parents: dict[str, list[Any]] = defaultdict(list)
    for edge in descriptor.relationships:
        parents[edge.child_member_id].append(edge)
    evidence = []
    for child_id, edges in sorted(parents.items()):
        ordered_edges = sorted(edges, key=lambda edge: (edge.parent_member_id, edge.assertion_id))
        members = [
            by_member[child_id],
            *(by_member[edge.parent_member_id] for edge in ordered_edges),
        ]
        variants: dict[str, dict[str, Any]] = {}
        for item in members:
            for row in item.variants:
                existing = variants.get(str(row["variant_id"]))
                canonical = {
                    key: row[key]
                    for key in ("assembly", "chromosome", "position", "reference", "alternate")
                }
                if existing is not None and any(
                    existing[key] != canonical[key] for key in canonical
                ):
                    raise ValueError("canonical variant identity conflicts across M2 runs")
                variants[str(row["variant_id"])] = {"variant_id": row["variant_id"], **canonical}
        genotype_maps: list[dict[str, list[dict[str, Any]]]] = []
        for item in members:
            grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
            for row in item.genotypes:
                grouped[str(row["variant_id"])].append(row)
            genotype_maps.append(grouped)
        for variant_id, variant in sorted(variants.items()):
            rows = [group.get(variant_id, []) for group in genotype_maps]
            conflicting = any(
                len(group) > 1 and len({tuple(row["alleles"]) for row in group}) > 1
                for group in rows
            )
            unsupported_ploidy = any(row.get("ploidy") != 2 for group in rows for row in group)

            def call(group: list[dict[str, Any]]) -> tuple[str, str] | None:
                distinct = {tuple(row["alleles"]) for row in group}
                if not group or len(distinct) != 1 or group[0].get("ploidy") != 2:
                    return None
                alleles = next(iter(distinct))
                return (str(alleles[0]), str(alleles[1]))

            evidence.append(
                evaluate_site(
                    variant_id=variant_id,
                    child_member_id=child_id,
                    child_genotype=call(rows[0]),
                    parent_genotypes=tuple(call(group) for group in rows[1:]),
                    relationship_assertion_ids=tuple(edge.assertion_id for edge in ordered_edges),
                    genotype_record_ids=tuple(
                        sorted(str(row["genotype_id"]) for group in rows for row in group)
                    ),
                    reference=str(variant["reference"]),
                    alternate=str(variant["alternate"]),
                    chromosome=str(variant["chromosome"]),
                    conflicting=conflicting,
                    unsupported_ploidy=unsupported_ploidy,
                )
            )
    counts = dict(sorted(Counter(item.compatibility.value for item in evidence).items()))
    stage = output.with_name(f".{output.name}.{run_id}.tmp")
    stage.mkdir(parents=True, mode=0o700)
    try:
        descriptor_copy = descriptor.model_dump(mode="json", by_alias=True)
        (stage / "pedigree.json").write_bytes(_dump(descriptor_copy))
        subject_rows = [
            {
                "member_id": item.member.member_id,
                "subject_id": item.member.subject_id,
                "m2_run_id": item.run_id,
                "m2_manifest_sha256": item.manifest_sha256,
            }
            for item in inputs
        ]
        _table(stage / "subject_inputs.parquet", subject_rows, {})
        evidence_rows = [
            {
                **item.model_dump(mode="json", exclude={"assignments"}),
                "assignment_count": len(item.assignments),
            }
            for item in evidence
        ]
        _table(stage / "segregation_evidence.parquet", evidence_rows, {})
        assignment_rows = [
            {
                "evidence_id": item.evidence_id,
                "transmission_status": item.transmission_status.value,
                "assignment_index": index,
                **assignment.model_dump(),
            }
            for item in evidence
            for index, assignment in enumerate(item.assignments)
        ]
        _table(
            stage / "transmission_phase_evidence.parquet",
            assignment_rows,
            {
                "evidence_id": pl.String,
                "transmission_status": pl.String,
                "assignment_index": pl.Int64,
                "parent_1_allele": pl.String,
                "parent_2_allele": pl.String,
            },
        )
        qc_rows = [
            {"code": status.value, "count": counts.get(status.value, 0)}
            for status in CompatibilityStatus
        ]
        _table(stage / "qc_findings.parquet", qc_rows, {})
        summary = {
            "schema": RUN_SCHEMA,
            "site_relationship_evidence_count": len(evidence),
            "compatibility_counts": counts,
            "compatibility_denominator": len(evidence),
            "rates_reported": False,
        }
        report = (
            "# Family segregation evidence\n\n"
            "> Conditional on declared relationships and directly observed M2 calls. "
            "This does not verify relatedness, call de novo variants, or provide "
            "clinical interpretation.\n\n"
            f"- Evaluated site/child relationship rows: {len(evidence)}\n"
            + "".join(f"- {key}: {value}\n" for key, value in counts.items())
        )
        (stage / "summary.json").write_bytes(_dump(summary))
        (stage / "report.md").write_text(report)
        metadata = {
            "schema": RUN_SCHEMA,
            "run_id": run_id,
            "package_version": __version__,
            "git_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], capture_output=True, text=True
            ).stdout.strip(),
            "created_at": datetime.now(UTC).isoformat(),
            "pedigree_source_file_sha256": _hash(pedigree),
            **scientific,
        }
        (stage / "family_analysis_metadata.json").write_bytes(_dump(metadata))
        artifact_names = sorted(path.name for path in stage.iterdir())
        artifacts = {
            name: {"sha256": _hash(stage / name), "byte_size": (stage / name).stat().st_size}
            for name in artifact_names
        }
        (stage / "manifest.json").write_bytes(
            _dump({"schema": RUN_SCHEMA, "run_id": run_id, "artifacts": artifacts})
        )
        (stage / "COMPLETED.json").write_bytes(
            _dump(
                {
                    "schema": "genome-evidence-m9-family-completion/v1",
                    "run_id": run_id,
                    "manifest_sha256": _hash(stage / "manifest.json"),
                }
            )
        )
        stage.rename(output)
        for path in output.iterdir():
            if path.is_file():
                os.chmod(path, 0o600)
        os.chmod(output, 0o700)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return FamilyAnalysisResult(
        run_id=run_id,
        output_directory=output,
        evidence_count=len(evidence),
        compatibility_counts=counts,
    )
