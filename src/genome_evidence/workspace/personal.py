"""Fail-closed personal notebook orchestration over the private workspace."""

import json
import os
import re
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

from genome_evidence.ingest import Ingest23andMeConfig, ingest_23andme
from genome_evidence.normalization import NormalizationConfig, normalize_m1_run
from genome_evidence.normalization.resources import canonical_assembly

from .core import (
    WorkspaceConfig,
    initialize_workspace,
    publish_completed_run,
    resolve_completed_run,
    resolve_latest_compatible_run,
    validate_workspace,
)

SOURCE_DIGEST = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class PersonalNormalizationResult:
    """Privacy-safe aggregate result from one published personal M1→M2 execution."""

    m1_run: Path
    m2_run: Path
    m1_run_id: str
    m2_run_id: str
    observation_count: int
    mapping_count: int
    canonical_genotype_count: int


def _hash_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _workspace_file(root: Path, value: str, label: str) -> Path:
    candidate = Path(value).expanduser()
    candidate = candidate if candidate.is_absolute() else root / candidate
    if candidate.is_symlink():
        raise ValueError(f"{label} must be a regular file inside the private workspace")
    candidate = candidate.resolve()
    if (
        not candidate.is_relative_to(root.resolve())
        or not candidate.is_file()
        or candidate.is_symlink()
    ):
        raise ValueError(f"{label} must be a regular file inside the private workspace")
    return candidate


def _single_resource(
    root: Path,
    environment: Mapping[str, str],
    variable: str,
    directory: str,
    suffixes: tuple[str, ...],
    label: str,
) -> Path:
    explicit = environment.get(variable)
    if explicit:
        return _workspace_file(root, explicit, label)
    resource_root = root / directory
    candidates = sorted(
        path.resolve()
        for path in resource_root.rglob("*")
        if path.is_file()
        and not path.is_symlink()
        and path.suffix.lower() in suffixes
        and not path.name.endswith(".fai")
    )
    if len(candidates) != 1:
        raise ValueError(
            f"{label} missing or ambiguous; install exactly one under {directory} "
            f"or set {variable} to a workspace-contained file"
        )
    return candidates[0]


def _selected_source(root: Path, environment: Mapping[str, str]) -> tuple[Path, str]:
    requested = environment.get("GENOME_EVIDENCE_SOURCE_SHA256")
    if requested is not None and not SOURCE_DIGEST.fullmatch(requested):
        raise ValueError("GENOME_EVIDENCE_SOURCE_SHA256 must be a lowercase SHA-256 digest")
    rows: list[tuple[Path, str]] = []
    workspace_config = WorkspaceConfig.model_validate_json(
        (root / "config/workspace.json").read_bytes()
    )
    for directory in sorted((root / "inputs/raw/23andme").iterdir()):
        if not directory.is_dir() or directory.is_symlink():
            continue
        try:
            manifest = json.loads((directory / "source_manifest.json").read_text())
        except (OSError, json.JSONDecodeError):
            continue
        digest = manifest.get("sha256")
        source = directory / "genome.txt"
        if (
            isinstance(digest, str)
            and digest == directory.name
            and SOURCE_DIGEST.fullmatch(digest)
            and manifest.get("data_file") == "genome.txt"
            and manifest.get("subject_id") == workspace_config.subject_id
            and source.is_file()
            and not source.is_symlink()
            and source.stat().st_size == manifest.get("byte_size")
            and _hash_file(source) == digest
        ):
            rows.append((source, digest))
    selected = [row for row in rows if requested is None or row[1] == requested]
    if len(selected) != 1:
        raise ValueError(
            "private source missing or ambiguous; complete notebook 00 and, when multiple "
            "content-addressed sources exist, set GENOME_EVIDENCE_SOURCE_SHA256"
        )
    return selected[0]


def resolve_personal_m2_run(root: Path, environment: Mapping[str, str] | None = None) -> Path:
    """Resolve an explicit or latest checksum-valid completed GRCh38 M2 run."""
    root = validate_workspace(root)
    env = os.environ if environment is None else environment
    constraints = {"target_assembly": "GRCh38"}
    explicit = env.get("GENOME_EVIDENCE_NORMALIZATION_RUN")
    if explicit:
        selected = Path(explicit).expanduser()
        selected = selected if selected.is_absolute() else root / selected
        return resolve_completed_run(root, selected, "M2", constraints)
    return resolve_latest_compatible_run(root, "M2", constraints)


def resolve_personal_population_bundle(
    root: Path, environment: Mapping[str, str] | None = None
) -> Path:
    """Resolve one explicitly reviewed local M5 bundle without downloading anything."""
    root = validate_workspace(root)
    env = os.environ if environment is None else environment
    explicit = env.get("GENOME_EVIDENCE_POPULATION_BUNDLE")
    if explicit:
        selected = Path(explicit).expanduser()
        selected = selected if selected.is_absolute() else root / selected
        if selected.is_symlink():
            raise ValueError("M5 reference bundle must be a directory inside the private workspace")
        selected = selected.resolve()
        if not selected.is_relative_to(root) or not selected.is_dir() or selected.is_symlink():
            raise ValueError("M5 reference bundle must be a directory inside the private workspace")
        return selected
    candidates = sorted(
        path.parent.resolve()
        for path in (root / "references/population_structure").rglob("manifest.json")
        if path.is_file() and not path.is_symlink() and not path.parent.is_symlink()
    )
    if len(candidates) != 1:
        raise ValueError(
            "reviewed M5 bundle missing or ambiguous; install exactly one bundle under "
            "references/population_structure or set GENOME_EVIDENCE_POPULATION_BUNDLE. "
            "The repository does not ship a production population reference bundle."
        )
    return candidates[0]


def run_personal_m1_m2(
    root: Path,
    subject_id: str,
    environment: Mapping[str, str] | None = None,
    *,
    working_root: Path | None = None,
) -> PersonalNormalizationResult:
    """Ingest one imported source, normalize it, and publish verified M1/M2 runs."""
    env = os.environ if environment is None else environment
    root = initialize_workspace(root, WorkspaceConfig(subject_id=subject_id))
    source, source_digest = _selected_source(root, env)
    markers = _single_resource(
        root,
        env,
        "GENOME_EVIDENCE_MARKER_DEFINITIONS",
        "references/markers/23andme",
        (".json",),
        "23andMe marker definitions",
    )
    fasta = _single_resource(
        root,
        env,
        "GENOME_EVIDENCE_GRCH38_FASTA",
        "references/genome/grch38",
        (".fa", ".fasta", ".fna"),
        "GRCh38 FASTA",
    )
    source_build_override = env.get("GENOME_EVIDENCE_SOURCE_BUILD")
    if working_root is not None:
        working_root = working_root.expanduser().resolve()
        working_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="genome-evidence-m1-m2-", dir=working_root) as temp:
        work = Path(temp)
        m1 = ingest_23andme(
            source,
            work / "m1",
            Ingest23andMeConfig(
                genome_build_override=source_build_override,
                sample_id=subject_id,
                source_identifier=f"workspace-sha256-{source_digest[:12]}",
            ),
        )
        resolved_source = canonical_assembly(m1.source_metadata.resolved_build)
        if resolved_source is None:
            raise ValueError(
                "source genome build is missing or unsupported; set "
                "GENOME_EVIDENCE_SOURCE_BUILD to GRCh37 or GRCh38 only after verification"
            )
        published_m1 = publish_completed_run(
            root,
            m1.output_directory,
            "M1",
            {"source_sha256": source_digest, "resolved_build": resolved_source},
        )
        liftover = None
        if resolved_source != "GRCh38":
            liftover = _single_resource(
                root,
                env,
                "GENOME_EVIDENCE_GRCH37_TO_GRCH38_LIFTOVER",
                "references/liftover/grch37_to_grch38",
                (".json",),
                "GRCh37-to-GRCh38 liftover map",
            )
        m2 = normalize_m1_run(
            m1.output_directory,
            work / "m2",
            NormalizationConfig(
                marker_definitions=markers,
                target_reference=fasta,
                marker_version=env.get("GENOME_EVIDENCE_MARKER_VERSION", "unversioned-local"),
                reference_version=env.get("GENOME_EVIDENCE_REFERENCE_VERSION", "GRCh38-local"),
                liftover=liftover,
                liftover_version=env.get("GENOME_EVIDENCE_LIFTOVER_VERSION", "unversioned-local"),
                source_build_override=source_build_override,
            ),
        )
        metadata = json.loads((m2.output_directory / "normalization_metadata.json").read_text())
        resource_hashes = {
            resource["resource_type"]: resource["sha256"] for resource in metadata["resources"]
        }
        constraints = {
            "source_sha256": source_digest,
            "target_assembly": "GRCh38",
            "marker_sha256": resource_hashes["marker_definitions"],
            "reference_sha256": resource_hashes["reference_sequence"],
        }
        if liftover is not None:
            constraints["liftover_sha256"] = resource_hashes["liftover"]
        published_m2 = publish_completed_run(root, m2.output_directory, "M2", constraints)
    return PersonalNormalizationResult(
        m1_run=published_m1,
        m2_run=published_m2,
        m1_run_id=m1.run_id,
        m2_run_id=m2.run_id,
        observation_count=len(m1.observations),
        mapping_count=len(m2.mappings),
        canonical_genotype_count=len(m2.genotypes),
    )
