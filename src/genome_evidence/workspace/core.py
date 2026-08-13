"""Fail-closed filesystem workspace with content-addressed private inputs."""

import json
import os
import re
import shutil
from hashlib import sha256
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

SCHEMA = "genome-evidence-workspace/v1"
SUBJECT = re.compile(r"^subject-[0-9]{4,}$")
DIRECTORIES = (
    "config/profiles",
    "inputs/inbox/23andme",
    "inputs/raw/23andme",
    "references/genome/grch38",
    "references/markers/23andme",
    "references/clinvar",
    "references/population_structure",
    "references/phasing_imputation",
    "references/polygenic_scores",
    "cache/downloads",
    "cache/tools/beagle",
    "cache/transformed_references",
    "cache/pgs_catalog",
    "registry/runs",
    "registry/latest",
    "runs/m1_ingestion",
    "runs/m2_normalization",
    "runs/m3_evidence",
    "runs/m3_linking",
    "runs/m4_prioritization",
    "runs/m5_population_structure",
    "runs/m6_phase_impute/_incomplete",
    "runs/m6_phase_impute/completed",
    "runs/m7_polygenic_scores/_incomplete",
    "runs/m7_polygenic_scores/completed",
    "reports",
    "exports",
    "logs",
)


class WorkspaceConfig(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", populate_by_name=True)
    subject_id: str = Field(pattern=r"^subject-[0-9]{4,}$")
    schema_id: str = Field(default=SCHEMA, alias="schema")


def _json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _inside_repo(root: Path) -> bool:
    resolved = root.resolve()
    return any((parent / ".git").exists() for parent in (resolved, *resolved.parents))


def initialize_workspace(root: Path, config: WorkspaceConfig) -> Path:
    """Create the canonical tree idempotently, never replacing conflicting config."""
    root = root.expanduser().resolve()
    if _inside_repo(root):
        raise ValueError("workspace must not be inside a Git checkout")
    root.mkdir(parents=True, exist_ok=True)
    for relative in DIRECTORIES:
        (root / relative).mkdir(parents=True, exist_ok=True)
    path = root / "config/workspace.json"
    data = config.model_dump(mode="json", by_alias=True)
    if path.exists():
        if json.loads(path.read_text(encoding="utf-8")) != data:
            raise FileExistsError("conflicting workspace configuration")
    else:
        _json(path, data)
    return validate_workspace(root)


def validate_workspace(root: Path) -> Path:
    root = root.expanduser().resolve()
    try:
        WorkspaceConfig.model_validate_json((root / "config/workspace.json").read_bytes())
    except Exception as error:
        raise ValueError("invalid or unknown workspace configuration") from error
    missing = [name for name in DIRECTORIES if not (root / name).is_dir()]
    if missing:
        raise ValueError(f"workspace directories missing: {', '.join(missing)}")
    return root


def import_23andme_source(root: Path, inbox_path: Path | None, subject_id: str) -> Path:
    """Copy exact bytes into immutable SHA-256 storage; return the source directory."""
    root = validate_workspace(root)
    config = WorkspaceConfig.model_validate_json((root / "config/workspace.json").read_bytes())
    if subject_id != config.subject_id or not SUBJECT.fullmatch(subject_id):
        raise ValueError("subject does not match the pseudonymous workspace subject")
    inbox = root / "inputs/inbox/23andme"
    if inbox_path is None:
        candidates = [p for p in inbox.iterdir() if p.is_file() and not p.is_symlink()]
        if len(candidates) != 1:
            raise ValueError(
                "inbox must contain exactly one regular file; select --file explicitly"
            )
        source = candidates[0]
    else:
        source = inbox_path.expanduser().resolve()
        if source.parent != inbox.resolve() or not source.is_file() or source.is_symlink():
            raise ValueError("source must be a regular file directly inside inputs/inbox/23andme")
    digest = sha256(source.read_bytes()).hexdigest()
    destination = root / "inputs/raw/23andme" / digest
    manifest = {
        "schema": "genome-evidence-private-source/v1",
        "subject_id": subject_id,
        "sha256": digest,
        "byte_size": source.stat().st_size,
        "data_file": "genome.txt",
        "privacy_class": "private_genotype_source",
    }
    if destination.exists():
        if sha256((destination / "genome.txt").read_bytes()).hexdigest() != digest:
            raise ValueError("immutable source destination conflicts with its identity")
        if json.loads((destination / "source_manifest.json").read_text()) != manifest:
            raise ValueError("immutable source manifest conflict")
        return destination
    destination.mkdir()
    shutil.copyfile(source, destination / "genome.txt")
    os.chmod(destination / "genome.txt", 0o600)
    _json(destination / "source_manifest.json", manifest)
    return destination


def _safe_relative(root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("registry paths must be workspace-relative")
    resolved = (root / path).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError("registry path escapes workspace")
    return resolved


def register_completed_run(root: Path, completion: dict[str, Any]) -> Path:
    """Register a verified immutable completion and update a JSON latest pointer."""
    root = validate_workspace(root)
    required = {"schema", "milestone", "run_id", "subject_id", "path", "constraints"}
    if set(completion) != required or completion["schema"] != "genome-evidence-completion/v1":
        raise ValueError("invalid completion schema")
    run = _safe_relative(root, str(completion["path"]))
    marker = run / "COMPLETED.json"
    if not marker.is_file() or json.loads(marker.read_text()) != completion:
        raise ValueError("completion marker missing or does not match registration")
    registry = root / "registry/runs" / f"{completion['run_id']}.json"
    if registry.exists() and json.loads(registry.read_text()) != completion:
        raise FileExistsError("immutable run registration conflict")
    if not registry.exists():
        _json(registry, completion)
    _json(
        root / "registry/latest" / f"{completion['milestone']}.json",
        {
            "schema": "genome-evidence-latest/v1",
            "run_id": completion["run_id"],
            "registry": str(registry.relative_to(root)),
        },
    )
    return registry


def list_completed_runs(root: Path, milestone: str | None = None) -> tuple[dict[str, Any], ...]:
    root = validate_workspace(root)
    rows = []
    for path in sorted((root / "registry/runs").glob("*.json")):
        row = json.loads(path.read_text())
        if milestone is None or row.get("milestone") == milestone:
            rows.append(row)
    return tuple(rows)


def resolve_latest_compatible_run(root: Path, milestone: str, constraints: dict[str, str]) -> Path:
    root = validate_workspace(root)
    pointer = json.loads((root / "registry/latest" / f"{milestone}.json").read_text())
    if pointer.get("schema") != "genome-evidence-latest/v1":
        raise ValueError("unknown latest pointer schema")
    registry = _safe_relative(root, pointer["registry"])
    completion = json.loads(registry.read_text())
    if completion.get("constraints") != constraints:
        raise ValueError("latest run is incompatible with requested constraints")
    run = _safe_relative(root, completion["path"])
    if json.loads((run / "COMPLETED.json").read_text()) != completion:
        raise ValueError("completed run failed registry validation")
    return run
