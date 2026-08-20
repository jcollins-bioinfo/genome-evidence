"""Checksum-valid M2 loading with explicit M9 subject binding.

M2 records predate a durable subject field, so M9 uses the narrow descriptor
binding and verifies it against the workspace completion marker.  Only canonical
genotypes from the M2 artifact are loaded; M6 records have no accepted interface.
"""

import json
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, cast

import polars as pl

from .models import PedigreeMember


@dataclass(frozen=True)
class M2Input:
    """Validated rows and identities for one explicitly bound subject."""

    member: PedigreeMember
    run_id: str
    manifest_sha256: str
    manifest: dict[str, Any]
    metadata: dict[str, Any]
    variants: tuple[dict[str, Any], ...]
    genotypes: tuple[dict[str, Any], ...]


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_m2_input(member: PedigreeMember, descriptor_directory: Path) -> M2Input:
    """Resolve and verify one completed M2 run without filename inference."""
    run = member.m2_run
    run = run if run.is_absolute() else descriptor_directory / run
    run = run.expanduser().resolve()
    if run.is_symlink() or not run.is_dir():
        raise ValueError("bound M2 run is missing or unsafe")
    try:
        completion = cast(dict[str, Any], json.loads((run / "COMPLETED.json").read_text()))
        manifest = cast(dict[str, Any], json.loads((run / "manifest.json").read_text()))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("bound M2 run lacks a valid completion marker or manifest") from error
    if (
        completion.get("schema") != "genome-evidence-completion/v1"
        or completion.get("milestone") != "M2"
        or completion.get("run_id") != manifest.get("run_id")
        or completion.get("subject_id") != member.subject_id
    ):
        raise ValueError("M2 completion identity does not match the explicit subject binding")
    if manifest.get("schema_version") != 1:
        raise ValueError("incompatible M2 schema")
    required = {"variants.parquet", "canonical_genotypes.parquet", "normalization_metadata.json"}
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or not required <= set(artifacts):
        raise ValueError("incomplete M2 artifact inventory")
    declared = {"manifest.json", "COMPLETED.json"}
    for name, identity in artifacts.items():
        if not isinstance(name, str) or Path(name).name != name:
            raise ValueError("unsafe M2 artifact name")
        path = run / name
        expected = identity if isinstance(identity, str) else identity.get("sha256")
        if not path.is_file() or path.is_symlink() or not isinstance(expected, str):
            raise ValueError("missing or invalid M2 artifact")
        if _hash(path) != expected:
            raise ValueError("M2 artifact checksum mismatch")
        declared.add(name)
    actual = {path.name for path in run.iterdir() if path.is_file()}
    if actual != declared or any(path.is_symlink() for path in run.iterdir()):
        raise ValueError("M2 run contains undeclared or unsafe files")
    metadata = cast(dict[str, Any], json.loads((run / "normalization_metadata.json").read_text()))
    if (
        metadata.get("run_id") != manifest.get("run_id")
        or metadata.get("target_assembly") != "GRCh38"
        or metadata.get("algorithm") != "m2-snv-1"
    ):
        raise ValueError("M2 identity, assembly, or canonical representation is incompatible")
    variants = tuple(pl.read_parquet(run / "variants.parquet").to_dicts())
    genotypes = tuple(pl.read_parquet(run / "canonical_genotypes.parquet").to_dicts())
    if any(row.get("normalization_run_id") != manifest.get("run_id") for row in genotypes):
        raise ValueError("M2 genotype lineage is incompatible with the selected run")
    return M2Input(
        member=member,
        run_id=str(manifest["run_id"]),
        manifest_sha256=_hash(run / "manifest.json"),
        manifest=manifest,
        metadata=metadata,
        variants=variants,
        genotypes=genotypes,
    )
