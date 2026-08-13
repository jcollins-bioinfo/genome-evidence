"""Validation for immutable, checksummed local PGS bundles."""

import json
from decimal import Decimal, InvalidOperation
from hashlib import sha256
from pathlib import Path
from typing import Any

import polars as pl

from .models import BundleValidation

BUNDLE_SCHEMA = "genome-evidence-pgs-bundle/v1"
REQUIRED = {"model_metadata.json", "model_variants.parquet", "model_exclusions.parquet"}


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_polygenic_score_bundle(bundle_directory: Path) -> BundleValidation:
    """Fail closed unless the complete local bundle and narrow model contract validate."""
    root = bundle_directory.resolve()
    manifest_path = root / "manifest.json"
    try:
        manifest: dict[str, Any] = json.loads(manifest_path.read_text())
    except Exception as error:
        raise ValueError("missing or invalid PGS bundle manifest.json") from error
    if manifest.get("schema") != BUNDLE_SCHEMA or set(manifest.get("artifacts", {})) != REQUIRED:
        raise ValueError("unknown bundle schema or incomplete artifact inventory")
    for relative, item in manifest["artifacts"].items():
        path = root / relative
        if Path(relative).is_absolute() or ".." in Path(relative).parts:
            raise ValueError("bundle artifact paths must be relative")
        if (
            not path.is_file()
            or path.stat().st_size != item.get("byte_size")
            or _hash(path) != item.get("sha256")
        ):
            raise ValueError(f"bundle artifact integrity failure: {relative}")
        if item.get("privacy_class") != "public_or_controlled_aggregate":
            raise ValueError("bundle must not contain target data")
    metadata = json.loads((root / "model_metadata.json").read_text())
    models = metadata.get("models", [])
    ids = tuple(x.get("pgs_id", "") for x in models)
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("model identities must be nonempty and unique")
    for model in models:
        required = {
            "pgs_id",
            "version",
            "trait",
            "source_url",
            "citation",
            "license",
            "assembly",
            "declared_variant_count",
            "completeness_policy",
        }
        if (
            not required <= set(model)
            or model["assembly"] != "GRCh38"
            or not str(model["source_url"]).startswith("https://")
        ):
            raise ValueError("incomplete or unsupported model provenance")
    rows = pl.read_parquet(root / "model_variants.parquet").to_dicts()
    keys: list[tuple[str, str, int, str, str, str]] = []
    for row in rows:
        key = (
            row["pgs_id"],
            row["assembly"],
            row["chromosome"],
            row["position"],
            row["reference"],
            row["alternate"],
        )
        if (
            row["pgs_id"] not in ids
            or row["assembly"] != "GRCh38"
            or row["chromosome"] not in {str(x) for x in range(1, 23)}
        ):
            raise ValueError("variant outside supported model/build/chromosome scope")
        if (
            len(row["reference"]) != 1
            or len(row["alternate"]) != 1
            or set(row["reference"] + row["alternate"]) - set("ACGT")
            or row["reference"] == row["alternate"]
        ):
            raise ValueError("only canonical biallelic A/C/G/T SNVs are supported")
        if row["effect_allele"] not in (row["reference"], row["alternate"]):
            raise ValueError("effect allele must equal canonical REF or ALT")
        try:
            weight = Decimal(row["effect_weight"])
        except (InvalidOperation, TypeError) as error:
            raise ValueError("invalid deterministic effect weight") from error
        if not weight.is_finite():
            raise ValueError("effect weight must be finite")
        keys.append(key)
    if keys != sorted(keys) or len(keys) != len(set(keys)):
        raise ValueError("model variants must be canonically sorted and unique")
    exclusions = pl.read_parquet(root / "model_exclusions.parquet")
    for model in models:
        prepared = sum(row["pgs_id"] == model["pgs_id"] for row in rows)
        excluded = (
            exclusions.filter(pl.col("pgs_id") == model["pgs_id"]).height
            if "pgs_id" in exclusions.columns
            else 0
        )
        if prepared + excluded != model["declared_variant_count"]:
            raise ValueError("declared/prepared/excluded count mismatch")
    return BundleValidation(
        directory=root, bundle_id=manifest["bundle_id"], model_ids=ids, artifact_count=len(REQUIRED)
    )
