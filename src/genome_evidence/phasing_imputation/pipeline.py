"""Offline M6 orchestration and strict M2 preflight."""

import json
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
from typing import Any

import polars as pl

from .engine import BeagleEngine
from .models import M6Config
from .reference import validate_phasing_reference


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_m2(root: Path) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    manifest_path = root / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("invalid M2 manifest") from error
    required = {
        "variants.parquet",
        "observation_mappings.parquet",
        "canonical_genotypes.parquet",
        "mapping_candidates.parquet",
        "normalization_qc.json",
        "normalization_report.md",
        "normalization_metadata.json",
    }
    if manifest.get("schema_version") != 1 or set(manifest.get("artifacts", {})) != required:
        raise ValueError("incomplete or unknown M2 artifact inventory")
    for name, digest in manifest["artifacts"].items():
        if not (root / name).is_file() or _hash(root / name) != digest:
            raise ValueError(f"M2 artifact integrity failure: {name}")
    metadata = json.loads((root / "normalization_metadata.json").read_text())
    if (
        metadata.get("run_id") != manifest.get("run_id")
        or metadata.get("target_assembly") != "GRCh38"
    ):
        raise ValueError("M2 run/assembly lineage mismatch")
    try:
        variants = pl.read_parquet(root / "variants.parquet").to_dicts()
        genotypes = pl.read_parquet(root / "canonical_genotypes.parquet").to_dicts()
    except Exception as error:
        raise ValueError("invalid M2 typed artifact schema") from error
    variant_ids = [str(x["variant_id"]) for x in variants]
    genotype_ids = [str(x["genotype_id"]) for x in genotypes]
    if len(set(variant_ids)) != len(variant_ids) or len(set(genotype_ids)) != len(genotype_ids):
        raise ValueError("M2 scientific identifiers are not unique")
    if any(str(x["variant_id"]) not in set(variant_ids) for x in genotypes):
        raise ValueError("M2 genotype has no canonical variant")
    return metadata, variants, genotypes


def align_m2(
    root: Path, chromosomes: tuple[str, ...]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Select only exact eligible calls; collapse only agreeing duplicate allele multisets."""
    _, variants, genotypes = _load_m2(root)
    by_variant = {str(v["variant_id"]): v for v in variants}
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for genotype in genotypes:
        grouped[str(genotype["variant_id"])].append(genotype)
    aligned, excluded = [], []
    bases = set("ACGT")
    for variant_id in sorted(grouped):
        variant = by_variant[variant_id]
        calls = grouped[variant_id]
        reason = None
        chromosome = str(variant["chromosome"])
        reference, alternate = str(variant["reference_allele"]), str(variant["alternate_allele"])
        allele_sets = {tuple(sorted(str(a) for a in call["alleles"])) for call in calls}
        if len(allele_sets) != 1:
            reason = "discordant_duplicate"
        elif chromosome not in chromosomes or chromosome not in {str(i) for i in range(1, 23)}:
            reason = "unsupported_chromosome"
        elif any(int(call["ploidy"]) != 2 for call in calls):
            reason = "unsupported_ploidy"
        elif len(reference) != 1 or len(alternate) != 1:
            reason = "non_snv"
        elif reference not in bases or alternate not in bases or reference == alternate:
            reason = "invalid_allele"
        elif any(set(call["alleles"]) - {reference, alternate} for call in calls):
            reason = "reference_mismatch"
        if reason:
            excluded.append(
                {
                    "variant_id": variant_id,
                    "reason": reason,
                    "genotype_ids": [x["genotype_id"] for x in calls],
                }
            )
        else:
            aligned.append(
                {
                    "variant_id": variant_id,
                    "assembly": "GRCh38",
                    "chromosome": chromosome,
                    "position": int(variant["position"]),
                    "reference": reference,
                    "alternate": alternate,
                    "alleles": list(next(iter(allele_sets))),
                    "genotype_ids": [x["genotype_id"] for x in calls],
                }
            )
    return sorted(
        aligned, key=lambda x: (int(x["chromosome"]), x["position"], x["reference"], x["alternate"])
    ), excluded


def phase_and_impute(
    normalization_directory: Path,
    reference_bundle_directory: Path,
    engine: BeagleEngine,
    output_directory: Path,
    config: M6Config,
) -> None:
    """Validate all inputs before an offline run; production execution is intentionally gated.

    The adapter and contracts are implemented in M6, but this release refuses to claim a
    completed scientific result until the reviewed production reference is installed and the
    real-backend smoke is enabled.
    """
    _load_m2(normalization_directory)
    reference = validate_phasing_reference(reference_bundle_directory)
    if not set(config.chromosomes).issubset(reference.manifest.chromosomes):
        raise ValueError("requested chromosome is absent from reference bundle")
    if output_directory.exists():
        raise FileExistsError("M6 never overwrites an output directory")
    align_m2(normalization_directory, config.chromosomes)
    if not engine.jar.is_file():
        raise ValueError("pinned local Beagle JAR is required; analysis never downloads it")
    raise RuntimeError(
        "M6 execution is explicitly incomplete: real Beagle fabricated-panel smoke has not "
        "been established; no COMPLETED.json was published"
    )
