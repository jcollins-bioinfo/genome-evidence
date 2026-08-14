"""Offline validation for immutable, checksummed M8 reference bundles.

Validation reads public/reference metadata only, rejects undeclared content and
symlinks, checks every byte, then enforces referential and scientific-scope
invariants.  It performs no network access and never accepts target genotype data.
"""

import json
from hashlib import sha256
from pathlib import Path
from typing import Any

import polars as pl

from .models import BundleValidation, ConstraintState, GeneStrategy

BUNDLE_SCHEMA = "genome-evidence-pgx-bundle/v1"
REQUIRED_ARTIFACTS = {
    "sources.json",
    "genes.parquet",
    "loci.parquet",
    "alleles.parquet",
    "allele_locus_constraints.parquet",
    "allele_function_assertions.parquet",
    "diplotype_phenotype_rules.parquet",
    "guideline_evidence.parquet",
}


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rows(root: Path, name: str) -> list[dict[str, Any]]:
    return pl.read_parquet(root / name).to_dicts()


def validate_pharmacogenomics_bundle(bundle_directory: Path) -> BundleValidation:
    """Validate a complete local M8 bundle without target data or network access.

    Parameters
    ----------
    bundle_directory:
        Directory containing ``bundle_manifest.json`` and its exact inventory.

    Returns
    -------
    BundleValidation
        Content-derived identity and supported gene inventory.

    Raises
    ------
    ValueError
        If schemas, paths, checksums, sources, ordering, foreign keys, constraints,
        capability declarations, or bounded candidate spaces are unsafe.

    Notes
    -----
    Unknown schema majors fail closed. Package and bundle schema versions are
    intentionally independent. Validation is deterministic and O(total rows).
    """
    root = bundle_directory.resolve()
    manifest_path = root / "bundle_manifest.json"
    try:
        manifest: dict[str, Any] = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as error:
        raise ValueError("missing or invalid pharmacogenomics bundle manifest") from error
    artifacts = manifest.get("artifacts", {})
    if manifest.get("schema") != BUNDLE_SCHEMA or set(artifacts) != REQUIRED_ARTIFACTS:
        raise ValueError("unknown bundle schema or incomplete artifact inventory")
    declared = {"bundle_manifest.json", *REQUIRED_ARTIFACTS}
    actual = {str(path.relative_to(root)) for path in root.iterdir() if path.name != "LICENSES"}
    if actual != declared:
        raise ValueError("bundle contains missing or extra top-level artifacts")
    for relative, identity in artifacts.items():
        rel = Path(relative)
        path = root / rel
        if rel.is_absolute() or ".." in rel.parts or path.is_symlink():
            raise ValueError("bundle artifact path is unsafe")
        if not path.is_file() or not path.resolve().is_relative_to(root):
            raise ValueError("bundle artifact is missing or escapes bundle")
        if path.stat().st_size != identity.get("byte_size") or _hash(path) != identity.get(
            "sha256"
        ):
            raise ValueError("bundle artifact checksum or size mismatch")
    sources = json.loads((root / "sources.json").read_text())
    source_rows = sources.get("sources", [])
    source_ids = {row.get("source_id") for row in source_rows}
    required_source = {
        "source_id",
        "resource",
        "version",
        "url",
        "retrieved_at",
        "sha256",
        "license",
        "content_fingerprint",
    }
    if (
        not source_rows
        or len(source_ids) != len(source_rows)
        or any(
            not required_source <= set(row)
            or not str(row["url"]).startswith("https://")
            or len(str(row["sha256"])) != 64
            or not row["version"]
            or not row["license"]
            for row in source_rows
        )
    ):
        raise ValueError(
            "source records require unique reproducible version/hash/license provenance"
        )
    genes = _rows(root, "genes.parquet")
    loci = _rows(root, "loci.parquet")
    alleles = _rows(root, "alleles.parquet")
    constraints = _rows(root, "allele_locus_constraints.parquet")
    functions = _rows(root, "allele_function_assertions.parquet")
    rules = _rows(root, "diplotype_phenotype_rules.parquet")
    guidelines = _rows(root, "guideline_evidence.parquet")
    for rows, id_column in (
        (genes, "gene_id"),
        (loci, "locus_id"),
        (alleles, "allele_id"),
    ):
        values = [str(row[id_column]) for row in rows]
        if not values or values != sorted(values) or len(values) != len(set(values)):
            raise ValueError(f"{id_column} values must be unique and canonically ordered")
    gene_ids = {row["gene_id"] for row in genes}
    locus_ids = {row["locus_id"] for row in loci}
    allele_ids = {row["allele_id"] for row in alleles}
    if any(row["gene_id"] not in gene_ids for row in loci + alleles + guidelines):
        raise ValueError("broken gene foreign key")
    if any(row["source_id"] not in source_ids for row in alleles + functions + rules + guidelines):
        raise ValueError("broken source foreign key")
    for locus in loci:
        if locus["assembly"] != "GRCh38" or locus["chromosome"] not in {
            str(i) for i in range(1, 23)
        }:
            raise ValueError("locus outside supported assembly/chromosome scope")
        ref, alt = locus["reference"], locus["alternate"]
        if not ref or not alt or ref == alt or set(ref + alt) - set("ACGT"):
            raise ValueError("invalid canonical REF/ALT representation")
    seen_constraints: set[tuple[str, str]] = set()
    for row in constraints:
        constraint_key = (row["allele_id"], row["locus_id"])
        if (
            constraint_key in seen_constraints
            or row["allele_id"] not in allele_ids
            or row["locus_id"] not in locus_ids
        ):
            raise ValueError("duplicate constraint or broken constraint foreign key")
        seen_constraints.add(constraint_key)
        ConstraintState(row["state"])
    if any(not any(a == row["allele_id"] for a, _ in seen_constraints) for row in alleles):
        raise ValueError("alleles cannot be defined through implicit absence")
    retired = {row["allele_id"] for row in alleles if row["status"] == "retired"}
    for row in rules:
        if row["allele_a"] not in allele_ids or row["allele_b"] not in allele_ids:
            raise ValueError("phenotype rule references undefined allele")
        if (row["allele_a"] in retired or row["allele_b"] in retired) and not row.get(
            "historical_mapping"
        ):
            raise ValueError("retired allele rule requires explicit historical mapping")
    strategies = {row["gene_id"]: GeneStrategy(row["strategy"]) for row in genes}
    structural = {
        row["allele_id"]
        for row in constraints
        if row["state"] == ConstraintState.UNSUPPORTED_STRUCTURAL_REQUIREMENT
    }
    if any(
        strategies[row["gene_id"]] == GeneStrategy.STAR_HAPLOTYPE_SMALL_VARIANT
        and row["allele_id"] in structural
        for row in alleles
    ):
        raise ValueError("generic matcher cannot ignore unsupported allele components")
    releases: dict[tuple[str, str], set[str]] = {}
    for row in alleles + functions + rules:
        releases.setdefault((row["gene_id"], row["source_id"]), set()).add(row["source_version"])
    if any(len(value) != 1 for value in releases.values()):
        raise ValueError("unsafely mixed source releases")
    digest = sha256(manifest_path.read_bytes()).hexdigest()
    return BundleValidation(
        directory=root,
        bundle_id=manifest["bundle_id"],
        bundle_version=manifest["bundle_version"],
        bundle_hash=digest,
        genes=tuple(sorted(gene_ids)),
    )
