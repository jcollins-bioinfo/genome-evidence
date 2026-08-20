"""Fail-closed filesystem workspace with content-addressed private inputs."""

import json
import os
import re
import shutil
from hashlib import sha256
from pathlib import Path
from typing import Any, cast

from pydantic import BaseModel, ConfigDict, Field

SCHEMA = "genome-evidence-workspace/v1"
SUBJECT = re.compile(r"^subject-[0-9]{4,}$")
DIRECTORIES = (
    "config/profiles",
    "inputs/inbox/23andme",
    "inputs/raw/23andme",
    "inputs/families",
    "references/genome/grch38",
    "references/markers/23andme",
    "references/liftover/grch37_to_grch38",
    "references/clinvar",
    "references/population_structure",
    "references/phasing_imputation",
    "references/polygenic_scores",
    "references/pharmacogenomics",
    "references/manifests/normalization",
    "cache/downloads",
    "cache/tools/beagle",
    "cache/tools/ucsc/kent-v479",
    "cache/transformed_references",
    "cache/pgs_catalog",
    "cache/clinpgx",
    "cache/pharmvar",
    "cache/pharmcat",
    "registry/runs",
    "registry/latest",
    "runs/m1_ingestion",
    "runs/m1_ingestion/_incomplete",
    "runs/m2_normalization",
    "runs/m2_normalization/_incomplete",
    "runs/m3_evidence",
    "runs/m3_linking",
    "runs/m4_prioritization",
    "runs/m5_population_structure",
    "runs/m5_population_structure/_incomplete",
    "runs/m6_phase_impute/_incomplete",
    "runs/m6_phase_impute/completed",
    "runs/m7_polygenic_scores/_incomplete",
    "runs/m7_polygenic_scores/completed",
    "runs/m8_pharmacogenomics/_incomplete",
    "runs/m8_pharmacogenomics/completed",
    "runs/m9_family_analysis/_incomplete",
    "runs/m9_family_analysis/completed",
    "reports",
    "exports",
    "logs",
)

RUN_DIRECTORIES = {
    "M1": "runs/m1_ingestion",
    "M2": "runs/m2_normalization",
    "M5": "runs/m5_population_structure",
}
RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


class WorkspaceConfig(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid", populate_by_name=True)
    subject_id: str = Field(pattern=r"^subject-[0-9]{4,}$")
    schema_id: str = Field(default=SCHEMA, alias="schema")


def _json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _inside_repo(root: Path) -> bool:
    resolved = root.resolve()
    for parent in (resolved, *resolved.parents):
        marker = parent / ".git"
        if marker.is_file() or (marker / "HEAD").is_file():
            return True
    return False


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


def _hash_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validated_run_manifest(run: Path, *, allow_completion: bool = False) -> dict[str, Any]:
    """Validate a run's declared files without reading scientific rows."""
    try:
        loaded = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("run manifest is missing or invalid") from error
    if not isinstance(loaded, dict):
        raise ValueError("run manifest must be a JSON object")
    manifest = cast(dict[str, Any], loaded)
    run_id = manifest.get("run_id")
    artifacts = manifest.get("artifacts")
    if not isinstance(run_id, str) or not RUN_ID.fullmatch(run_id):
        raise ValueError("run manifest has an invalid run ID")
    if not isinstance(artifacts, dict) or not artifacts:
        raise ValueError("run manifest has no artifact inventory")
    declared = {"manifest.json"}
    if allow_completion:
        declared.add("COMPLETED.json")
    for name, identity in artifacts.items():
        if not isinstance(name, str):
            raise ValueError("run artifact names must be strings")
        artifact = _safe_relative(run, name)
        if artifact.parent != run.resolve() or not artifact.is_file() or artifact.is_symlink():
            raise ValueError(f"run artifact is missing or unsafe: {name}")
        expected_hash: Any
        expected_size: Any
        if isinstance(identity, str):
            expected_hash, expected_size = identity, None
        elif isinstance(identity, dict):
            expected_hash = identity.get("sha256")
            expected_size = identity.get("byte_size")
        else:
            raise ValueError(f"run artifact identity is invalid: {name}")
        if not isinstance(expected_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            raise ValueError(f"run artifact checksum is invalid: {name}")
        if expected_size is not None and (
            not isinstance(expected_size, int) or artifact.stat().st_size != expected_size
        ):
            raise ValueError(f"run artifact size mismatch: {name}")
        if _hash_file(artifact) != expected_hash:
            raise ValueError(f"run artifact checksum mismatch: {name}")
        declared.add(name)
    actual = {
        str(path.relative_to(run))
        for path in run.rglob("*")
        if path.is_file() and not path.is_symlink()
    }
    if any(path.is_symlink() for path in run.rglob("*")):
        raise ValueError("run directories must not contain symlinks")
    if actual != declared:
        raise ValueError("run directory contains undeclared files")
    return manifest


def _validated_completion(
    root: Path, completion: dict[str, Any], *, expected_milestone: str | None = None
) -> Path:
    required = {"schema", "milestone", "run_id", "subject_id", "path", "constraints"}
    if set(completion) != required or completion.get("schema") != "genome-evidence-completion/v1":
        raise ValueError("invalid completion schema")
    milestone = completion.get("milestone")
    run_id = completion.get("run_id")
    constraints = completion.get("constraints")
    if milestone not in RUN_DIRECTORIES or (expected_milestone and milestone != expected_milestone):
        raise ValueError("unknown or incompatible completion milestone")
    if not isinstance(run_id, str) or not RUN_ID.fullmatch(run_id):
        raise ValueError("invalid completion run ID")
    if not isinstance(constraints, dict) or any(
        not isinstance(key, str) or not isinstance(value, str) for key, value in constraints.items()
    ):
        raise ValueError("completion constraints must be string pairs")
    config = WorkspaceConfig.model_validate_json((root / "config/workspace.json").read_bytes())
    if completion.get("subject_id") != config.subject_id:
        raise ValueError("completion subject does not match workspace")
    expected = (root / RUN_DIRECTORIES[milestone] / run_id).resolve()
    run = _safe_relative(root, str(completion.get("path")))
    if run != expected:
        raise ValueError("completion path does not match milestone and run ID")
    return run


def register_completed_run(root: Path, completion: dict[str, Any]) -> Path:
    """Register a verified immutable completion and update a JSON latest pointer."""
    root = validate_workspace(root)
    run = _validated_completion(root, completion)
    marker = run / "COMPLETED.json"
    if not marker.is_file() or json.loads(marker.read_text()) != completion:
        raise ValueError("completion marker missing or does not match registration")
    manifest = _validated_run_manifest(run, allow_completion=True)
    if manifest["run_id"] != completion["run_id"]:
        raise ValueError("completion and run manifest identities disagree")
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


def publish_completed_run(
    root: Path, source: Path, milestone: str, constraints: dict[str, str]
) -> Path:
    """Copy, re-hash, complete, and register one immutable private run."""
    root = validate_workspace(root)
    source = source.expanduser().resolve()
    if not source.is_dir() or source.is_symlink():
        raise ValueError("run source must be a regular directory")
    manifest = _validated_run_manifest(source)
    run_id = str(manifest["run_id"])
    if milestone not in RUN_DIRECTORIES:
        raise ValueError("unsupported workspace milestone")
    if any(
        not isinstance(key, str) or not isinstance(value, str) for key, value in constraints.items()
    ):
        raise ValueError("completion constraints must be string pairs")
    config = WorkspaceConfig.model_validate_json((root / "config/workspace.json").read_bytes())
    parent = root / RUN_DIRECTORIES[milestone]
    staging = parent / "_incomplete" / run_id
    destination = parent / run_id
    if staging.exists() or destination.exists():
        raise FileExistsError("workspace run destination already exists")
    completion = {
        "schema": "genome-evidence-completion/v1",
        "milestone": milestone,
        "run_id": run_id,
        "subject_id": config.subject_id,
        "path": str(destination.relative_to(root)),
        "constraints": constraints,
    }
    shutil.copytree(source, staging)
    try:
        _validated_run_manifest(staging)
        shutil.copytree(staging, destination)
        _validated_run_manifest(destination)
        _json(destination / "COMPLETED.json", completion)
        register_completed_run(root, completion)
    except Exception:
        if destination.exists() and not (destination / "COMPLETED.json").exists():
            shutil.rmtree(destination)
        raise
    else:
        shutil.rmtree(staging)
    return destination


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
    try:
        pointer = json.loads(
            (root / "registry/latest" / f"{milestone}.json").read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"no completed {milestone} run is registered") from error
    if pointer.get("schema") != "genome-evidence-latest/v1":
        raise ValueError("unknown latest pointer schema")
    registry = _safe_relative(root, pointer["registry"])
    completion = json.loads(registry.read_text())
    run = _validated_completion(root, completion, expected_milestone=milestone)
    if pointer.get("run_id") != completion.get("run_id"):
        raise ValueError("latest pointer and completion identities disagree")
    actual_constraints = completion["constraints"]
    if any(actual_constraints.get(key) != value for key, value in constraints.items()):
        raise ValueError("latest run is incompatible with requested constraints")
    if json.loads((run / "COMPLETED.json").read_text()) != completion:
        raise ValueError("completed run failed registry validation")
    _validated_run_manifest(run, allow_completion=True)
    return run


def resolve_completed_run(
    root: Path, run: Path, milestone: str, constraints: dict[str, str]
) -> Path:
    """Validate an explicitly selected completed workspace run."""
    root = validate_workspace(root)
    selected = run.expanduser().resolve()
    if not selected.is_relative_to(root):
        raise ValueError("selected run must be inside the private workspace")
    try:
        completion = json.loads((selected / "COMPLETED.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("selected run has no valid completion marker") from error
    validated = _validated_completion(root, completion, expected_milestone=milestone)
    if validated != selected:
        raise ValueError("selected run does not match its completion path")
    actual_constraints = completion["constraints"]
    if any(actual_constraints.get(key) != value for key, value in constraints.items()):
        raise ValueError("selected run is incompatible with requested constraints")
    _validated_run_manifest(selected, allow_completion=True)
    return selected
