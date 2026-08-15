"""Pinned, provenance-bearing normalization-resource provisioning for personal data."""

import gzip
import json
import os
import random
import re
import shutil
import subprocess
import tempfile
import threading
import time
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import UTC, datetime
from hashlib import md5, sha256
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from genome_evidence.normalization.resources import canonical_assembly, canonical_chromosome
from genome_evidence.workspace.provisioning_progress import (
    ProvisioningReporter,
    Sleeper,
    resumable_download,
)
from genome_evidence.workspace.segmented_download import segmented_download

DBSNP_BUILD = "155"
DBSNP_URLS = {
    "GRCh37": "https://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp155.bb",
    "GRCh38": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb",
}
DBSNP_QUERY_URLS = {assembly: (url,) for assembly, url in DBSNP_URLS.items()}
DBSNP_COMMON_URLS = {
    "GRCh37": "https://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp155Common.bb",
    "GRCh38": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb",
}
FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
FASTA_UPSTREAM_MD5 = "1c9dcaddfa41027f17cd8f7a82c7293b"
KENT_VERSION = "479"
KENT_TOOL_URL = (
    f"https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v{KENT_VERSION}/bigBedNamedItems"
)
_KENT_USAGE_TOKENS = (
    "bigbednameditems - extract item of given name from bigbed",
    "bigbednameditems file.bb name output.bed",
    "-namefile",
)
SELECTION_SCHEMA = "genome-evidence-normalization-resource-selection/v1"
PROVENANCE_SCHEMA = "genome-evidence-normalization-resource-provenance/v1"
CHECKPOINT_SCHEMA = "genome-evidence-normalization-resource-checkpoint/v1"
QUERY_CHECKPOINT_SCHEMA = "genome-evidence-dbsnp-query-checkpoint/v1"
QUERY_SPLIT_SCHEMA = "genome-evidence-dbsnp-query-split/v1"
BUNDLE_COMPLETION_SCHEMA = "genome-evidence-normalization-resource-bundle-completion/v1"
SELECTION_PUBLICATION_SCHEMA = "genome-evidence-normalization-selection-publication/v1"
SELECTION_COMPLETION_SCHEMA = "genome-evidence-normalization-selection-completion/v1"
RESOURCE_ALGORITHM_VERSION = "genome-evidence-normalization-resource-builder/v3-common-parallel"
_BUILD_PATTERN = re.compile(r"(?:build|assembly)[\s:=]+(GRCh\d+|hg\d+|37|38)\b", re.IGNORECASE)
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_RSID = re.compile(r"^rs\d+$")
_DEFAULT_BATCH_SIZE = 5_000
_DEFAULT_QUERY_ATTEMPTS = 4
_DEFAULT_DBSNP_WORKERS = 6
_DEFAULT_DOWNLOAD_SEGMENTS = 8
_MIN_SPLIT_BATCH_SIZE = 250


@dataclass(frozen=True)
class SourceMarker:
    marker_id: str
    chromosome: str
    position: int


@dataclass(frozen=True)
class BigBedVariant:
    chromosome: str
    position: int
    marker_id: str
    reference: str
    alternates: tuple[str, ...]


class NormalizationResourceSelection(BaseModel):
    """Durable, source-compatible pointers to one validated resource set."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    schema_id: str = Field(default=SELECTION_SCHEMA, alias="schema")
    source_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    source_assembly: str
    marker_definitions: str
    grch38_fasta: str
    grch37_to_grch38_liftover: str | None = None
    marker_version: str
    reference_version: str
    liftover_version: str | None = None
    provenance_manifest: str


@dataclass(frozen=True)
class ProvisioningResult:
    selection: NormalizationResourceSelection
    selection_path: Path
    marker_count: int
    mapped_marker_count: int
    unresolved_marker_count: int
    log_path: Path
    checkpoint_path: Path


class ProvisioningIncomplete(RuntimeError):
    """Operational interruption whose durable checkpoints can be resumed safely."""

    def __init__(self, message: str, *, log_path: Path, checkpoint_path: Path) -> None:
        super().__init__(message)
        self.log_path = log_path
        self.checkpoint_path = checkpoint_path


class _OperationalInterruption(RuntimeError):
    """Internal marker for retry-exhausted operations that remain safely resumable."""


def _selection_control_paths(root: Path) -> tuple[Path, Path, Path]:
    selection = root / "config/normalization_resources.json"
    publishing = root / "config/normalization_resources.PUBLISHING.json"
    completed = root / "config/normalization_resources.COMPLETED.json"
    return selection, publishing, completed


def _selection_completion_matches(root: Path, source_sha256: str) -> bool:
    selection, _, completed = _selection_control_paths(root)
    if completed.is_symlink():
        return False
    try:
        manifest = json.loads(completed.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    return bool(
        manifest.get("schema") == SELECTION_COMPLETION_SCHEMA
        and manifest.get("source_sha256") == source_sha256
        and selection.is_file()
        and not selection.is_symlink()
        and selection.stat().st_size == manifest.get("byte_size")
        and _hash(selection) == manifest.get("sha256")
    )


def load_normalization_resource_selection(
    root: Path,
    source_sha256: str,
    *,
    reporter: ProvisioningReporter | None = None,
) -> NormalizationResourceSelection | None:
    """Load a durable selection and prove that it belongs to the selected private source."""
    path, _, completed_path = _selection_control_paths(root)
    if not path.exists():
        if completed_path.exists():
            raise ValueError("normalization resource selection completion marker is orphaned")
        return None
    if not path.is_file() or path.is_symlink():
        raise ValueError("normalization resource selection path is unsafe")
    if completed_path.exists() and not _selection_completion_matches(root, source_sha256):
        raise ValueError("normalization resource selection completion marker is invalid")
    try:
        selection = NormalizationResourceSelection.model_validate_json(path.read_bytes())
    except Exception as error:
        raise ValueError("normalization resource selection is invalid") from error
    if selection.source_sha256 != source_sha256:
        raise ValueError(
            "normalization resources were provisioned for a different source; rerun notebook 00B"
        )
    if canonical_assembly(selection.source_assembly) != selection.source_assembly:
        raise ValueError("normalization resource selection has an unsupported source assembly")
    required_artifacts = {
        selection.marker_definitions,
        selection.grch38_fasta,
        f"{selection.grch38_fasta}.fai",
    }
    if selection.grch37_to_grch38_liftover is not None:
        required_artifacts.add(selection.grch37_to_grch38_liftover)
    provenance_relative = Path(selection.provenance_manifest)
    if provenance_relative.is_absolute() or ".." in provenance_relative.parts:
        raise ValueError("normalization resource selection paths must be workspace-relative")
    provenance_path = (root / provenance_relative).resolve()
    if (
        not provenance_path.is_relative_to(root.resolve())
        or not provenance_path.is_file()
        or provenance_path.is_symlink()
    ):
        raise ValueError("normalization resource selection references a missing or unsafe file")
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("normalization resource provenance manifest is invalid") from error
    if (
        provenance.get("schema") != PROVENANCE_SCHEMA
        or provenance.get("source_sha256") != source_sha256
        or provenance.get("source_assembly") != selection.source_assembly
    ):
        raise ValueError("normalization resource provenance is incompatible with the selection")
    artifacts = provenance.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ValueError("normalization resource provenance has no artifact inventory")
    if not required_artifacts.issubset(artifacts):
        raise ValueError("normalization resource provenance does not cover every selection")
    if completed_path.exists():
        selection_completion = json.loads(completed_path.read_text(encoding="utf-8"))
        bundle_id = provenance.get("bundle_id")
        if not isinstance(bundle_id, str) or selection_completion.get("bundle_id") != bundle_id:
            raise ValueError("normalization resource selector names an incompatible bundle")
        bundle_completion_path = provenance_path.with_name(f"{bundle_id}.COMPLETED.json")
        bundle_completion = _load_bundle_completion(
            root,
            bundle_completion_path,
            bundle_id=bundle_id,
            source_sha256=source_sha256,
        )
        if bundle_completion is None:
            raise ValueError("normalization resource bundle completion marker is invalid")
        if (
            bundle_completion.get("artifacts") != artifacts
            or bundle_completion.get("provenance", {}).get("path") != selection.provenance_manifest
        ):
            raise ValueError("normalization resource bundle completion inventory conflicts")
        return selection
    for index, (relative, identity) in enumerate(sorted(artifacts.items()), start=1):
        if not isinstance(identity, dict) or not _DIGEST.fullmatch(str(identity.get("sha256"))):
            raise ValueError("normalization resource provenance has an invalid artifact identity")
        candidate = Path(relative)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError("normalization resource provenance paths must be workspace-relative")
        resource = (root / candidate).resolve()
        if (
            not resource.is_relative_to(root.resolve())
            or not resource.is_file()
            or resource.is_symlink()
        ):
            raise ValueError(
                "normalization resource provenance references a missing or unsafe file"
            )
        if resource.stat().st_size != identity.get("byte_size"):
            raise ValueError("normalization resource size does not match provenance")
        if (
            _hash(
                resource,
                reporter=reporter,
                event=f"selection.verify.{index}",
            )
            != identity["sha256"]
        ):
            raise ValueError("normalization resource checksum does not match provenance")
    return selection


Download = Callable[[str, Path], None]
Runner = Callable[..., subprocess.CompletedProcess[str]]


def _hash(
    path: Path,
    algorithm: str = "sha256",
    *,
    reporter: ProvisioningReporter | None = None,
    event: str | None = None,
) -> str:
    digest = sha256() if algorithm == "sha256" else md5()  # noqa: S324 - upstream publishes MD5
    started_at = time.monotonic()
    total = path.stat().st_size
    completed = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            completed += len(chunk)
            if reporter is not None and event is not None:
                reporter.progress(
                    event,
                    f"Hashing {path.name} ({algorithm.upper()})",
                    completed,
                    total,
                    unit="bytes",
                    started_at=started_at,
                    algorithm=algorithm,
                )
    if reporter is not None and event is not None:
        reporter.progress(
            event,
            f"Hashing {path.name} ({algorithm.upper()})",
            completed,
            total,
            unit="bytes",
            started_at=started_at,
            force=True,
            algorithm=algorithm,
        )
    return digest.hexdigest()


def _obtain_download(
    url: str,
    destination: Path,
    download: Download | None,
    reporter: ProvisioningReporter,
    *,
    sleep: Sleeper,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if download is None:
        try:
            resumable_download(url, destination, reporter, sleep=sleep)
        except (OSError, RuntimeError) as error:
            raise _OperationalInterruption(str(error)) from error
        return
    reporter.info(
        "download.injected",
        f"Acquiring {destination.name} through the configured transfer adapter.",
        url=url,
        path=str(destination),
    )
    download(url, destination)
    if not destination.is_file():
        raise RuntimeError("configured transfer adapter did not create the requested file")
    reporter.success(
        "download.complete",
        f"Acquired {destination.name} ({destination.stat().st_size:,} bytes)",
        url=url,
        path=str(destination),
        byte_size=destination.stat().st_size,
    )


def _json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(payload, encoding="utf-8")
    os.replace(temporary, path)


def _workspace_relative(root: Path, path: Path) -> str:
    resolved = path.resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError("resource path escapes the private workspace")
    return str(resolved.relative_to(root.resolve()))


def _load_bundle_completion(
    root: Path,
    path: Path,
    *,
    bundle_id: str,
    source_sha256: str,
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    if path.is_symlink():
        raise ValueError("normalization resource bundle completion marker is unsafe")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    if not isinstance(manifest, dict):
        return None
    if (
        manifest.get("schema") != BUNDLE_COMPLETION_SCHEMA
        or manifest.get("bundle_id") != bundle_id
        or manifest.get("source_sha256") != source_sha256
        or manifest.get("algorithm_version") != RESOURCE_ALGORITHM_VERSION
    ):
        raise ValueError("normalization resource bundle completion identity is invalid")
    identities = manifest.get("artifacts")
    provenance = manifest.get("provenance")
    if not isinstance(identities, dict) or not isinstance(provenance, dict):
        raise ValueError("normalization resource bundle completion inventory is invalid")
    provenance_path = provenance.get("path")
    if not isinstance(provenance_path, str):
        raise ValueError("normalization resource bundle completion provenance is invalid")
    inventory = {**identities, provenance_path: provenance}
    for relative, identity in inventory.items():
        if not isinstance(identity, dict) or not _DIGEST.fullmatch(str(identity.get("sha256"))):
            raise ValueError("normalization resource bundle completion checksum is invalid")
        candidate = Path(relative)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError("normalization resource bundle completion path is unsafe")
        artifact = (root / candidate).resolve()
        if (
            not artifact.is_relative_to(root.resolve())
            or not artifact.is_file()
            or artifact.is_symlink()
            or artifact.stat().st_size != identity.get("byte_size")
            or _hash(artifact) != identity.get("sha256")
        ):
            raise ValueError("completed normalization resource bundle failed validation")
    return manifest


def _selected_source(root: Path, environment: Mapping[str, str]) -> tuple[Path, str]:
    requested = environment.get("GENOME_EVIDENCE_SOURCE_SHA256")
    if requested is not None and not _DIGEST.fullmatch(requested):
        raise ValueError("GENOME_EVIDENCE_SOURCE_SHA256 must be a lowercase SHA-256 digest")
    candidates: list[tuple[Path, str]] = []
    for directory in sorted((root / "inputs/raw/23andme").iterdir()):
        if (
            not directory.is_dir()
            or directory.is_symlink()
            or not _DIGEST.fullmatch(directory.name)
        ):
            continue
        source = directory / "genome.txt"
        manifest_path = directory / "source_manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        if (
            source.is_file()
            and not source.is_symlink()
            and manifest.get("sha256") == directory.name
            and manifest.get("data_file") == "genome.txt"
            and source.stat().st_size == manifest.get("byte_size")
            and _hash(source) == directory.name
        ):
            candidates.append((source, directory.name))
    selected = [row for row in candidates if requested is None or row[1] == requested]
    if len(selected) != 1:
        raise ValueError(
            "private source missing or ambiguous; complete notebook 00 and select one source digest"
        )
    return selected[0]


def read_source_markers(
    path: Path, build_override: str | None = None
) -> tuple[str, tuple[SourceMarker, ...]]:
    """Read identifiers and coordinates only; genotype tokens are neither retained nor emitted."""
    try:
        text = path.read_bytes().decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ValueError("23andMe source must be UTF-8 text") from error
    comments: list[str] = []
    markers: list[SourceMarker] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        if line.startswith("#"):
            comments.append(line)
            continue
        columns = line.split("\t")
        if len(columns) != 4:
            raise ValueError(f"invalid source record at line {line_number}")
        marker_id, chromosome, position_text, _genotype = columns
        try:
            position = int(position_text)
        except ValueError as error:
            raise ValueError(f"invalid source coordinate at line {line_number}") from error
        if position <= 0:
            raise ValueError(f"invalid source coordinate at line {line_number}")
        markers.append(SourceMarker(marker_id, canonical_chromosome(chromosome), position))
    declared = next(
        (match.group(1) for comment in comments if (match := _BUILD_PATTERN.search(comment))), None
    )
    resolved = canonical_assembly(build_override or declared or "")
    if resolved not in {"GRCh37", "GRCh38"}:
        raise ValueError(
            "source build is missing or unsupported; set GENOME_EVIDENCE_SOURCE_BUILD only "
            "after independently verifying GRCh37 or GRCh38"
        )
    if not markers:
        raise ValueError("23andMe source contains no marker records")
    return resolved, tuple(markers)


def parse_bigbed_variants(text: str) -> tuple[BigBedVariant, ...]:
    """Parse the documented bigDbSnp BED4+ fields used for exact rsID extraction."""
    rows: list[BigBedVariant] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line:
            continue
        columns = line.split("\t")
        if len(columns) < 7:
            raise ValueError(f"invalid bigDbSnp row at line {line_number}")
        try:
            start, end, alt_count = int(columns[1]), int(columns[2]), int(columns[5])
        except ValueError as error:
            raise ValueError(f"invalid bigDbSnp coordinates at line {line_number}") from error
        reference = columns[4].upper()
        alternates = tuple(x.upper() for x in columns[6].rstrip(",").split(",") if x)
        if end <= start or alt_count != len(alternates):
            raise ValueError(f"incoherent bigDbSnp row at line {line_number}")
        rows.append(
            BigBedVariant(
                canonical_chromosome(columns[0]),
                start + 1,
                columns[3],
                reference,
                alternates,
            )
        )
    return tuple(rows)


def build_marker_resources(
    markers: Sequence[SourceMarker],
    source_assembly: str,
    source_variants: Sequence[BigBedVariant],
    target_variants: Sequence[BigBedVariant],
) -> tuple[list[dict[str, Any]], dict[str, list[list[str | int]]], dict[str, int]]:
    """Create exact-coordinate SNV definitions and conservative cross-build mappings."""
    marker_index = {
        (x.marker_id, x.chromosome, x.position) for x in markers if x.marker_id.startswith("rs")
    }
    accepted_source = [
        row
        for row in source_variants
        if (row.marker_id, row.chromosome, row.position) in marker_index
    ]
    target_by_id: dict[str, list[BigBedVariant]] = {}
    for row in target_variants:
        target_by_id.setdefault(row.marker_id, []).append(row)
    definitions: list[dict[str, Any]] = []
    definition_keys: set[tuple[str, str, int, str, str]] = set()
    mapping_sets: dict[str, set[tuple[str, int]]] = {}
    mapped_ids: set[str] = set()
    for row in sorted(
        accepted_source, key=lambda x: (x.chromosome, x.position, x.marker_id, x.reference)
    ):
        valid_alts = sorted(
            {
                alt
                for alt in row.alternates
                if len(row.reference) == len(alt) == 1
                and set(row.reference + alt) <= set("ACGT")
                and alt != row.reference
            }
        )
        for alt in valid_alts:
            definition_key = (
                row.marker_id,
                row.chromosome,
                row.position,
                row.reference,
                alt,
            )
            if definition_key in definition_keys:
                continue
            definition_keys.add(definition_key)
            definitions.append(
                {
                    "marker_id": row.marker_id,
                    "assembly": source_assembly,
                    "chromosome": row.chromosome,
                    "position": row.position,
                    "reference": row.reference,
                    "alternate": alt,
                    "orientation": "none",
                    "orientation_authoritative": True,
                    "rsid": row.marker_id,
                }
            )
        if source_assembly == "GRCh37" and valid_alts:
            compatible = {
                (target.chromosome, target.position)
                for target in target_by_id.get(row.marker_id, [])
                if target.reference == row.reference
                and any(alt in target.alternates for alt in valid_alts)
                and len(target.reference) == 1
            }
            if compatible:
                key = f"{row.chromosome}:{row.position - 1}"
                mapping_sets.setdefault(key, set()).update(compatible)
                mapped_ids.add(row.marker_id)
    mappings: dict[str, list[list[str | int]]] = {
        key: [[chromosome, position - 1] for chromosome, position in sorted(values)]
        for key, values in sorted(mapping_sets.items())
    }
    source_ids = {row.marker_id for row in accepted_source}
    stats = {
        "source_marker_count": len(markers),
        "rsid_marker_count": len({x.marker_id for x in markers if x.marker_id.startswith("rs")}),
        "defined_marker_count": len({row["marker_id"] for row in definitions}),
        "definition_count": len(definitions),
        "cross_build_mapped_marker_count": len(mapped_ids),
        "exact_source_placement_count": len(source_ids),
    }
    return definitions, mappings, stats


def build_fasta_index(
    fasta: Path,
    output: Path,
    *,
    reporter: ProvisioningReporter | None = None,
    event: str = "fasta.index",
) -> None:
    """Build a standard .fai while rejecting irregular non-terminal sequence lines."""
    entries: list[tuple[str, int, int, int, int]] = []
    name: str | None = None
    length = 0
    offset = 0
    line_bases = 0
    line_width = 0
    previous_bases: int | None = None
    started_at = time.monotonic()
    total_bytes = fasta.stat().st_size
    with fasta.open("rb") as handle:
        while raw := handle.readline():
            if raw.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, line_bases, line_width))
                header = raw[1:].strip().split(maxsplit=1)[0]
                if not header:
                    raise ValueError("FASTA contains an empty sequence name")
                try:
                    name = header.decode("ascii")
                except UnicodeDecodeError as error:
                    raise ValueError("FASTA sequence names must be ASCII") from error
                length = 0
                offset = handle.tell()
                line_bases = line_width = 0
                previous_bases = None
                continue
            if name is None:
                raise ValueError("FASTA sequence data precedes the first header")
            bases = raw.rstrip(b"\r\n")
            if not bases:
                raise ValueError("FASTA contains a blank sequence line")
            if previous_bases is not None and previous_bases != line_bases:
                raise ValueError("FASTA has an irregular non-terminal sequence line")
            if line_bases == 0:
                line_bases, line_width = len(bases), len(raw)
            elif len(raw) != line_width and len(bases) == line_bases:
                raise ValueError("FASTA line terminators are inconsistent")
            previous_bases = len(bases)
            length += len(bases)
            if reporter is not None:
                reporter.progress(
                    event,
                    "Building hg38.fa.fai",
                    handle.tell(),
                    total_bytes,
                    unit="bytes",
                    started_at=started_at,
                )
    if name is not None:
        entries.append((name, length, offset, line_bases, line_width))
    if not entries or any(entry[1] == 0 for entry in entries):
        raise ValueError("FASTA contains no complete sequences")
    names = [entry[0] for entry in entries]
    if len(names) != len(set(names)):
        raise ValueError("FASTA contains duplicate sequence names")
    output.write_text(
        "".join("\t".join(map(str, entry)) + "\n" for entry in entries), encoding="ascii"
    )
    if reporter is not None:
        reporter.progress(
            event,
            "Building hg38.fa.fai",
            total_bytes,
            total_bytes,
            unit="bytes",
            started_at=started_at,
            force=True,
            sequence_count=len(entries),
        )


def _write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(value, encoding="ascii")
    os.replace(temporary, path)


def _safe_subprocess_diagnostic(value: str, identifier_path: Path) -> str:
    sanitized = value.replace(str(identifier_path), "[identifier-file]")
    sanitized = re.sub(r"\brs[\w.-]*", "[identifier]", sanitized, flags=re.IGNORECASE)
    return " ".join(sanitized.split())[-2000:]


def _valid_query_checkpoint(
    directory: Path,
    *,
    identifiers: Sequence[str],
    identifiers_sha256: str,
    assembly: str,
    canonical_url: str,
    tool_sha256: str,
) -> Path | None:
    output = directory / "records.bed"
    manifest_path = directory / "COMPLETED.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    if (
        manifest.get("schema") != QUERY_CHECKPOINT_SCHEMA
        or manifest.get("identifiers_sha256") != identifiers_sha256
        or manifest.get("identifier_count") != len(identifiers)
        or manifest.get("assembly") != assembly
        or manifest.get("dbsnp_build") != DBSNP_BUILD
        or manifest.get("canonical_url") != canonical_url
        or manifest.get("tool_sha256") != tool_sha256
        or not output.is_file()
        or output.is_symlink()
        or output.stat().st_size != manifest.get("byte_size")
        or _hash(output) != manifest.get("sha256")
    ):
        return None
    try:
        rows = parse_bigbed_variants(output.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError):
        return None
    requested = set(identifiers)
    if manifest.get("record_count") != len(rows) or any(
        row.marker_id not in requested for row in rows
    ):
        return None
    return output


def _publish_checkpoint_file(source: Path, destination: Path) -> None:
    """Copy verified local work into Drive; a completion sidecar is written separately."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as input_handle, destination.open("wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=4 * 1024 * 1024)
        output_handle.flush()
        os.fsync(output_handle.fileno())
    if source.stat().st_size != destination.stat().st_size or _hash(source) != _hash(destination):
        raise OSError(f"checkpoint copy verification failed for {destination.name}")


def _run_query_with_heartbeat(
    arguments: list[str],
    runner: Runner,
    reporter: ProvisioningReporter,
    *,
    identifier_count: int,
    attempt: int,
    split_depth: int,
    timeout_seconds: int,
    worker_tmp: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    process_environment = None
    if worker_tmp is not None:
        process_environment = {**os.environ, "TMPDIR": str(worker_tmp)}
    if runner is not subprocess.run:
        return runner(
            arguments,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            env=process_environment,
        )
    started_at = time.monotonic()
    process = subprocess.Popen(  # noqa: S603 - pinned local utility and fixed arguments
        arguments,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=process_environment,
    )
    try:
        while True:
            elapsed = time.monotonic() - started_at
            remaining = timeout_seconds - elapsed
            if remaining <= 0:
                process.kill()
                stdout, stderr = process.communicate()
                raise subprocess.TimeoutExpired(
                    arguments,
                    timeout_seconds,
                    output=stdout,
                    stderr=stderr,
                )
            try:
                stdout, stderr = process.communicate(timeout=min(30.0, remaining))
            except subprocess.TimeoutExpired:
                elapsed = time.monotonic() - started_at
                reporter.info(
                    "dbsnp.batch.heartbeat",
                    f"UCSC query is still active after {elapsed / 60:.1f} minutes.",
                    identifier_count=identifier_count,
                    attempt=attempt,
                    split_depth=split_depth,
                    elapsed_seconds=elapsed,
                    timeout_seconds=timeout_seconds,
                )
                continue
            return subprocess.CompletedProcess(
                arguments,
                process.returncode,
                stdout=stdout,
                stderr=stderr,
            )
    except BaseException:
        if process.poll() is None:
            process.kill()
            process.communicate()
        raise


def _query_batch(
    tool: Path,
    tool_sha256: str,
    urls: Sequence[str],
    assembly: str,
    identifiers: Sequence[str],
    directory: Path,
    work: Path,
    runner: Runner,
    reporter: ProvisioningReporter,
    *,
    attempts: int,
    timeout_seconds: int,
    sleep: Sleeper,
    split_depth: int = 0,
    worker_tmp: Path | None = None,
) -> Path:
    payload = "\n".join(identifiers) + "\n"
    identifiers_sha256 = sha256(payload.encode("ascii")).hexdigest()
    completed = _valid_query_checkpoint(
        directory,
        identifiers=identifiers,
        identifiers_sha256=identifiers_sha256,
        assembly=assembly,
        canonical_url=urls[0],
        tool_sha256=tool_sha256,
    )
    if completed is not None:
        reporter.info(
            "dbsnp.batch.resume",
            f"Reused a verified {len(identifiers):,}-identifier query checkpoint.",
            identifier_count=len(identifiers),
            checkpoint=str(directory),
        )
        return completed
    directory.mkdir(parents=True, exist_ok=True)
    scratch = work / "dbsnp-query" / identifiers_sha256
    scratch.mkdir(parents=True, exist_ok=True)
    identifier_path = scratch / "identifiers.txt"
    _write_text(identifier_path, payload)
    os.chmod(identifier_path, 0o600)
    output = directory / "records.bed"
    temporary = scratch / "records.bed.part"
    last_error = "no subprocess result"
    split_path = directory / "SPLIT.json"
    split_planned = False
    if len(identifiers) > _MIN_SPLIT_BATCH_SIZE and split_path.is_file():
        try:
            split_manifest = json.loads(split_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            split_manifest = {}
        split_planned = (
            split_manifest.get("schema") == QUERY_SPLIT_SCHEMA
            and split_manifest.get("assembly") == assembly
            and split_manifest.get("dbsnp_build") == DBSNP_BUILD
            and split_manifest.get("canonical_url") == urls[0]
            and split_manifest.get("tool_sha256") == tool_sha256
            and split_manifest.get("identifiers_sha256") == identifiers_sha256
            and split_manifest.get("identifier_count") == len(identifiers)
            and split_manifest.get("midpoint") == len(identifiers) // 2
            and split_manifest.get("minimum_batch_size") == _MIN_SPLIT_BATCH_SIZE
        )
        if split_planned:
            reporter.info(
                "dbsnp.batch.split.resume",
                "Resuming an adaptive split without repeating the failed parent query.",
                identifier_count=len(identifiers),
                split_depth=split_depth,
            )
    attempt_numbers = range(1, attempts + 1) if not split_planned else range(0)
    for attempt in attempt_numbers:
        temporary.unlink(missing_ok=True)
        endpoint = urls[(attempt - 1) % len(urls)]
        reporter.info(
            "dbsnp.batch.start",
            f"Querying {len(identifiers):,} identifiers (attempt {attempt}/{attempts}).",
            identifier_count=len(identifiers),
            attempt=attempt,
            attempts=attempts,
            endpoint=endpoint,
            split_depth=split_depth,
        )
        started_at = time.monotonic()
        try:
            process = _run_query_with_heartbeat(
                [str(tool), "-nameFile", endpoint, str(identifier_path), str(temporary)],
                runner,
                reporter,
                identifier_count=len(identifiers),
                attempt=attempt,
                split_depth=split_depth,
                timeout_seconds=timeout_seconds,
                worker_tmp=worker_tmp,
            )
            if process.returncode == 0 and temporary.is_file():
                validation_error_type: str | None = None
                try:
                    rows = parse_bigbed_variants(temporary.read_text(encoding="utf-8"))
                    requested = set(identifiers)
                    if any(row.marker_id not in requested for row in rows):
                        raise ValueError(
                            "UCSC query returned an identifier outside the requested batch"
                        )
                except (OSError, UnicodeError, ValueError) as error:
                    validation_error_type = type(error).__name__
                    last_error = (
                        "zero-exit query output failed validation "
                        f"({validation_error_type}); partial output was rejected"
                    )
                    rows = None
                if rows is None:
                    reporter.warning(
                        "dbsnp.batch.output-rejected",
                        "UCSC returned an invalid batch output; it will be discarded and retried.",
                        identifier_count=len(identifiers),
                        attempt=attempt,
                        split_depth=split_depth,
                        error_type=validation_error_type,
                    )
                    temporary.unlink(missing_ok=True)
                    if attempt < attempts:
                        delay = min(2 ** (attempt - 1), 30) + random.uniform(0.0, 1.0)
                        reporter.info(
                            "dbsnp.batch.backoff",
                            f"Retrying this query batch in {delay} seconds.",
                            identifier_count=len(identifiers),
                            attempt=attempt,
                            next_attempt=attempt + 1,
                            delay_seconds=delay,
                            split_depth=split_depth,
                        )
                        sleep(delay)
                    continue
                _publish_checkpoint_file(temporary, output)
                manifest = {
                    "schema": QUERY_CHECKPOINT_SCHEMA,
                    "assembly": assembly,
                    "dbsnp_build": DBSNP_BUILD,
                    "canonical_url": urls[0],
                    "endpoint_used": endpoint,
                    "tool_sha256": tool_sha256,
                    "identifiers_sha256": identifiers_sha256,
                    "identifier_count": len(identifiers),
                    "record_count": len(rows),
                    "sha256": _hash(output),
                    "byte_size": output.stat().st_size,
                }
                _json(directory / "COMPLETED.json", manifest)
                elapsed = max(time.monotonic() - started_at, 1e-9)
                reporter.success(
                    "dbsnp.batch.complete",
                    f"Query batch completed: {len(identifiers):,} identifiers, "
                    f"{len(rows):,} records, "
                    f"{len(identifiers) / elapsed:,.0f} identifiers/s.",
                    identifier_count=len(identifiers),
                    record_count=len(rows),
                    elapsed_seconds=elapsed,
                    identifiers_per_second=len(identifiers) / elapsed,
                    endpoint=endpoint,
                    split_depth=split_depth,
                )
                return output
            stderr = ((process.stderr or "") + (process.stdout or "")).strip()
            last_error = _safe_subprocess_diagnostic(stderr, identifier_path) or (
                f"exit status {process.returncode} without diagnostics"
            )
        except subprocess.TimeoutExpired as error:
            last_error = f"query timed out after {error.timeout} seconds"
        except subprocess.CalledProcessError as error:
            details = ((error.stderr or "") + (error.stdout or "")).strip()
            last_error = _safe_subprocess_diagnostic(details, identifier_path) or (
                f"exit status {error.returncode} without diagnostics"
            )
        except OSError as error:
            last_error = f"operating-system interruption: {type(error).__name__}"
        reporter.warning(
            "dbsnp.batch.retry",
            f"Query attempt failed; checkpoint remains resumable. Diagnostic: {last_error}",
            identifier_count=len(identifiers),
            attempt=attempt,
            attempts=attempts,
            endpoint=endpoint,
            split_depth=split_depth,
        )
        if attempt < attempts:
            delay = min(2 ** (attempt - 1), 30) + random.uniform(0.0, 1.0)
            reporter.info(
                "dbsnp.batch.backoff",
                f"Retrying this query batch in {delay} seconds.",
                identifier_count=len(identifiers),
                attempt=attempt,
                next_attempt=attempt + 1,
                delay_seconds=delay,
                split_depth=split_depth,
            )
            sleep(delay)
    if len(identifiers) > _MIN_SPLIT_BATCH_SIZE:
        midpoint = len(identifiers) // 2
        if not split_planned:
            _json(
                split_path,
                {
                    "schema": QUERY_SPLIT_SCHEMA,
                    "assembly": assembly,
                    "dbsnp_build": DBSNP_BUILD,
                    "canonical_url": urls[0],
                    "tool_sha256": tool_sha256,
                    "identifiers_sha256": identifiers_sha256,
                    "identifier_count": len(identifiers),
                    "midpoint": midpoint,
                    "minimum_batch_size": _MIN_SPLIT_BATCH_SIZE,
                },
            )
            reporter.warning(
                "dbsnp.batch.split",
                f"Repeated failure; bisecting {len(identifiers):,} identifiers into "
                "smaller resumable batches.",
                identifier_count=len(identifiers),
                split_depth=split_depth,
            )
        children = (
            _query_batch(
                tool,
                tool_sha256,
                urls,
                assembly,
                identifiers[:midpoint],
                directory / "split-left",
                work,
                runner,
                reporter,
                attempts=attempts,
                timeout_seconds=timeout_seconds,
                sleep=sleep,
                split_depth=split_depth + 1,
                worker_tmp=worker_tmp,
            ),
            _query_batch(
                tool,
                tool_sha256,
                urls,
                assembly,
                identifiers[midpoint:],
                directory / "split-right",
                work,
                runner,
                reporter,
                attempts=attempts,
                timeout_seconds=timeout_seconds,
                sleep=sleep,
                split_depth=split_depth + 1,
                worker_tmp=worker_tmp,
            ),
        )
        split_temporary = scratch / "split-records.bed.part"
        with split_temporary.open("wb") as destination:
            for child in children:
                with child.open("rb") as source:
                    shutil.copyfileobj(source, destination, length=1024 * 1024)
            destination.flush()
            os.fsync(destination.fileno())
        rows = parse_bigbed_variants(split_temporary.read_text(encoding="utf-8"))
        requested = set(identifiers)
        if any(row.marker_id not in requested for row in rows):
            raise ValueError("adaptive UCSC query returned an unrequested identifier")
        _publish_checkpoint_file(split_temporary, output)
        _json(
            directory / "COMPLETED.json",
            {
                "schema": QUERY_CHECKPOINT_SCHEMA,
                "assembly": assembly,
                "dbsnp_build": DBSNP_BUILD,
                "canonical_url": urls[0],
                "endpoint_used": "adaptive-split",
                "tool_sha256": tool_sha256,
                "identifiers_sha256": identifiers_sha256,
                "identifier_count": len(identifiers),
                "record_count": len(rows),
                "sha256": _hash(output),
                "byte_size": output.stat().st_size,
            },
        )
        return output
    raise _OperationalInterruption(
        f"UCSC dbSNP query failed after retries; rerun to resume from completed batches. "
        f"Last diagnostic: {last_error}"
    )


def _query_bigbed_in_batches(
    tool: Path,
    tool_sha256: str,
    urls: Sequence[str],
    assembly: str,
    identifiers: Sequence[str],
    checkpoint_root: Path,
    work: Path,
    runner: Runner,
    reporter: ProvisioningReporter,
    *,
    batch_size: int,
    attempts: int,
    timeout_seconds: int,
    sleep: Sleeper,
    label: str,
    workers: int = 1,
) -> Path:
    all_payload = "\n".join(identifiers) + "\n"
    all_sha256 = sha256(all_payload.encode("ascii")).hexdigest()
    query_key = sha256(
        f"{urls[0]}|{tool_sha256}|{batch_size}|{all_sha256}".encode("ascii")
    ).hexdigest()[:24]
    directory = checkpoint_root / label / query_key
    combined = _valid_query_checkpoint(
        directory,
        identifiers=identifiers,
        identifiers_sha256=all_sha256,
        assembly=assembly,
        canonical_url=urls[0],
        tool_sha256=tool_sha256,
    )
    if combined is not None:
        reporter.success(
            f"dbsnp.{label}.resume",
            f"Reused the complete verified {label} dbSNP extract "
            f"({len(identifiers):,} identifiers).",
            identifier_count=len(identifiers),
            checkpoint=str(directory),
        )
        return combined
    batches = [
        identifiers[index : index + batch_size] for index in range(0, len(identifiers), batch_size)
    ]
    reporter.info(
        f"dbsnp.{label}.plan",
        f"Planned {len(batches):,} resumable batches of at most {batch_size:,} identifiers.",
        identifier_count=len(identifiers),
        batch_count=len(batches),
        batch_size=batch_size,
        canonical_url=urls[0],
        checkpoint=str(directory),
    )
    batch_specs: list[tuple[int, Sequence[str], Path, Path | None]] = []
    for index, batch in enumerate(batches):
        batch_payload = "\n".join(batch) + "\n"
        batch_key = sha256(batch_payload.encode("ascii")).hexdigest()[:12]
        batch_directory = directory / "batches" / f"{index:05d}-{batch_key}"
        valid = _valid_query_checkpoint(
            batch_directory,
            identifiers=batch,
            identifiers_sha256=sha256(batch_payload.encode("ascii")).hexdigest(),
            assembly=assembly,
            canonical_url=urls[0],
            tool_sha256=tool_sha256,
        )
        batch_specs.append((index, batch, batch_directory, valid))
    resumed_identifiers = sum(len(batch) for _, batch, _, valid in batch_specs if valid)
    resumed_batches = sum(1 for *_, valid in batch_specs if valid)
    if resumed_batches:
        reporter.info(
            f"dbsnp.{label}.resume",
            f"Verified and resumed {resumed_batches:,}/{len(batches):,} completed batches.",
            resumed_batches=resumed_batches,
            resumed_identifiers=resumed_identifiers,
        )
    started_at = time.monotonic()
    completed_identifiers = resumed_identifiers
    completed_batches = resumed_batches
    outputs: list[Path | None] = [None] * len(batches)
    local = threading.local()
    worker_counter = iter(range(workers))
    worker_lock = threading.Lock()

    def run(spec: tuple[int, Sequence[str], Path, Path | None]) -> tuple[int, Path, int]:
        index, batch, batch_directory, valid = spec
        if valid is not None:
            return index, valid, 0
        if not hasattr(local, "work"):
            with worker_lock:
                worker_index = next(worker_counter)
            local.work = work / "kent-workers" / f"worker-{worker_index:02d}"
            local.work.mkdir(parents=True, exist_ok=True)
            (local.work / "udcCache").mkdir(exist_ok=True)
        result = _query_batch(
            tool,
            tool_sha256,
            urls,
            assembly,
            batch,
            batch_directory,
            local.work,
            runner,
            reporter,
            attempts=attempts,
            timeout_seconds=timeout_seconds,
            sleep=sleep,
            worker_tmp=local.work,
        )
        return index, result, len(batch)

    pending = [spec for spec in batch_specs if spec[3] is None]
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="dbsnp-fallback") as pool:
        futures = [pool.submit(run, spec) for spec in pending]
        try:
            for future in as_completed(futures):
                index, batch_output, count = future.result()
                outputs[index] = batch_output
                completed_identifiers += count
                completed_batches += 1
                reporter.progress(
                    f"dbsnp.{label}.progress",
                    f"dbSNP {label} query ({completed_batches:,}/{len(batches):,} batches)",
                    completed_identifiers,
                    len(identifiers),
                    unit="identifiers",
                    started_at=started_at,
                    initial_completed=resumed_identifiers,
                    force=True,
                    completed_batches=completed_batches,
                    total_batches=len(batches),
                    configured_workers=workers,
                    active_workers=min(workers, len(futures)),
                )
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    directory.mkdir(parents=True, exist_ok=True)
    temporary = work / "dbsnp-combine" / f"{query_key}.bed.part"
    temporary.parent.mkdir(parents=True, exist_ok=True)
    output = directory / "records.bed"
    with temporary.open("wb") as destination:
        for completed_output in outputs:
            if completed_output is None:
                raise RuntimeError("internal query batch output is missing")
            with completed_output.open("rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
        destination.flush()
        os.fsync(destination.fileno())
    rows = parse_bigbed_variants(temporary.read_text(encoding="utf-8"))
    requested = set(identifiers)
    if any(row.marker_id not in requested for row in rows):
        raise ValueError("combined UCSC query returned an unrequested identifier")
    _publish_checkpoint_file(temporary, output)
    _json(
        directory / "COMPLETED.json",
        {
            "schema": QUERY_CHECKPOINT_SCHEMA,
            "assembly": assembly,
            "dbsnp_build": DBSNP_BUILD,
            "canonical_url": urls[0],
            "endpoint_used": "batched-primary",
            "tool_sha256": tool_sha256,
            "identifiers_sha256": all_sha256,
            "identifier_count": len(identifiers),
            "batch_size": batch_size,
            "batch_count": len(batches),
            "record_count": len(rows),
            "sha256": _hash(output),
            "byte_size": output.stat().st_size,
        },
    )
    reporter.success(
        f"dbsnp.{label}.complete",
        f"Completed verified {label} extract with {len(rows):,} returned records.",
        identifier_count=len(identifiers),
        record_count=len(rows),
        batch_count=len(batches),
        output=str(output),
    )
    return output


def _merge_query_outputs(paths: Sequence[Path], destination: Path) -> Path:
    """Merge validated BED rows independent of worker completion order."""
    rows: set[str] = set()
    for path in paths:
        text = path.read_text(encoding="utf-8")
        parse_bigbed_variants(text)
        rows.update(line for line in text.splitlines() if line)
    destination.parent.mkdir(parents=True, exist_ok=True)
    _write_text(destination, "".join(f"{line}\n" for line in sorted(rows)))
    return destination


def _query_common_first(
    tool: Path,
    tool_sha256: str,
    assembly: str,
    identifiers: Sequence[str],
    checkpoint_root: Path,
    work: Path,
    runner: Runner,
    reporter: ProvisioningReporter,
    *,
    download: Download | None,
    batch_size: int,
    attempts: int,
    timeout_seconds: int,
    workers: int,
    download_segments: int,
    sleep: Sleeper,
    label: str,
) -> tuple[Path, dict[str, Any]]:
    """Use a local common BigBed as a non-authoritative accelerator, then full fallback."""
    common_cache = checkpoint_root / "common" / assembly.lower() / "dbSnp155Common.bb"
    common_manifest: dict[str, Any] | None = None
    common_output: Path | None = None
    try:
        if download is None:
            common_cache, common_manifest = segmented_download(
                DBSNP_COMMON_URLS[assembly],
                common_cache,
                reporter,
                concurrency=download_segments,
                attempts=attempts,
                timeout_seconds=timeout_seconds,
                sleep=sleep,
            )
        else:
            _obtain_download(
                DBSNP_COMMON_URLS[assembly], common_cache, download, reporter, sleep=sleep
            )
            common_manifest = {
                "remote_identity": {"canonical_url": DBSNP_COMMON_URLS[assembly]},
                "sha256": _hash(common_cache),
                "byte_size": common_cache.stat().st_size,
                "segment_concurrency": 0,
            }
        local_common = work / "common" / assembly.lower() / common_cache.name
        local_common.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(common_cache, local_common)
        if _hash(local_common) != common_manifest["sha256"]:
            raise OSError("local common dbSNP copy failed checksum verification")
        common_output = _query_bigbed_in_batches(
            tool,
            tool_sha256,
            (str(local_common),),
            assembly,
            identifiers,
            checkpoint_root / "common-query",
            work,
            runner,
            reporter,
            batch_size=batch_size,
            attempts=attempts,
            timeout_seconds=timeout_seconds,
            sleep=sleep,
            label=f"{label}-common",
            workers=1,
        )
        common_rows = parse_bigbed_variants(common_output.read_text(encoding="utf-8"))
        returned = {row.marker_id for row in common_rows}
        unresolved = tuple(identifier for identifier in identifiers if identifier not in returned)
        reporter.success(
            f"dbsnp.{label}.common.complete",
            f"Local common dbSNP returned {len(returned):,} identifiers; "
            f"{len(unresolved):,} require full fallback.",
            common_hit_count=len(returned),
            full_fallback_count=len(unresolved),
        )
    except Exception as error:
        reporter.warning(
            f"dbsnp.{label}.common.unavailable",
            "Common dbSNP could not be validated; safely querying every identifier "
            "against the full pinned index.",
            error_type=type(error).__name__,
            full_fallback_count=len(identifiers),
        )
        common_output = None
        unresolved = tuple(identifiers)
        common_manifest = {"status": "unavailable", "canonical_url": DBSNP_COMMON_URLS[assembly]}
    # v1 placed full-query batches directly below dbsnp/<assembly>/<query-key>. Rebuild
    # that deterministic plan and accept only independently valid completion manifests.
    legacy_outputs: list[Path] = []
    legacy_covered: set[str] = set()
    legacy_parent = checkpoint_root / label
    legacy_query_key = sha256(
        f"{DBSNP_URLS[assembly]}|{tool_sha256}|{batch_size}|"
        f"{sha256(('\n'.join(identifiers) + '\n').encode('ascii')).hexdigest()}".encode("ascii")
    ).hexdigest()[:24]
    for index, batch_start in enumerate(range(0, len(identifiers), batch_size)):
        batch = identifiers[batch_start : batch_start + batch_size]
        payload = "\n".join(batch) + "\n"
        batch_digest = sha256(payload.encode("ascii")).hexdigest()
        candidate = (
            legacy_parent / legacy_query_key / "batches" / f"{index:05d}-{batch_digest[:12]}"
        )
        valid = _valid_query_checkpoint(
            candidate,
            identifiers=batch,
            identifiers_sha256=batch_digest,
            assembly=assembly,
            canonical_url=DBSNP_URLS[assembly],
            tool_sha256=tool_sha256,
        )
        if valid is not None:
            legacy_outputs.append(valid)
            legacy_covered.update(batch)
    if legacy_covered:
        unresolved = tuple(
            identifier for identifier in unresolved if identifier not in legacy_covered
        )
        reporter.info(
            f"dbsnp.{label}.legacy.resume",
            f"Reused {len(legacy_outputs):,} verified legacy full-query batches.",
            resumed_legacy_batches=len(legacy_outputs),
            resumed_legacy_identifiers=len(legacy_covered),
        )
    fallback = None
    if unresolved:
        fallback = _query_bigbed_in_batches(
            tool,
            tool_sha256,
            DBSNP_QUERY_URLS[assembly],
            assembly,
            unresolved,
            checkpoint_root / "fallback",
            work,
            runner,
            reporter,
            batch_size=batch_size,
            attempts=attempts,
            timeout_seconds=timeout_seconds,
            sleep=sleep,
            label=f"{label}-full",
            workers=workers,
        )
        reporter.info(
            f"dbsnp.{label}.resume",
            "Verified fallback checkpoints were incorporated where available.",
            full_fallback_count=len(unresolved),
        )
    merged = checkpoint_root / "merged" / f"{label}.bed"
    paths = [path for path in (common_output, *legacy_outputs, fallback) if path is not None]
    _merge_query_outputs(paths, merged)
    return merged, {
        "common": common_manifest,
        "common_output": (
            {"sha256": _hash(common_output), "byte_size": common_output.stat().st_size}
            if common_output
            else None
        ),
        "full_url": DBSNP_URLS[assembly],
        "common_hit_count": len(identifiers) - len(unresolved),
        "full_fallback_count": len(unresolved),
        "legacy_batch_count": len(legacy_outputs),
        "workers": workers,
    }


def _prepare_kent_tool(
    root: Path,
    work: Path,
    download: Download | None,
    runner: Runner,
    reporter: ProvisioningReporter,
    *,
    sleep: Sleeper,
) -> tuple[Path, str]:
    cached = root / f"cache/tools/ucsc/kent-v{KENT_VERSION}/bigBedNamedItems"
    if not cached.is_file():
        cached.parent.mkdir(parents=True, exist_ok=True)
        _obtain_download(KENT_TOOL_URL, cached, download, reporter, sleep=sleep)
    if cached.is_symlink():
        raise ValueError("cached UCSC utility must not be a symlink")
    local = work / "bigBedNamedItems"
    for validation_attempt in range(1, 3):
        reporter.info(
            "kent.cache",
            f"Using pinned Kent v{KENT_VERSION} utility cache ({cached.stat().st_size:,} bytes).",
            path=str(cached),
            byte_size=cached.stat().st_size,
            pinned_version=KENT_VERSION,
        )
        shutil.copyfile(cached, local)
        os.chmod(local, 0o700)
        try:
            probe = runner([str(local)], capture_output=True, text=True, check=False, timeout=30)
        except OSError as error:
            probe = subprocess.CompletedProcess(
                [str(local)],
                126,
                stdout="",
                stderr=type(error).__name__,
            )
        probe_text = (probe.stdout or "") + (probe.stderr or "")
        normalized_probe = probe_text.lower()
        valid = probe.returncode != 0 and all(
            token in normalized_probe for token in _KENT_USAGE_TOKENS
        )
        if valid:
            break
        if download is not None or validation_attempt == 2:
            raise ValueError(
                "cached UCSC utility does not expose the expected bigBedNamedItems CLI"
            )
        invalid = cached.with_name(f"{cached.name}.invalid-{_hash(cached)[:12]}-{int(time.time())}")
        reporter.warning(
            "kent.cache.repair",
            "The cached Kent utility failed its CLI probe; preserving it as invalid and "
            "performing one clean download.",
            invalid_name=invalid.name,
        )
        os.replace(cached, invalid)
        cached.with_name(cached.name + ".part").unlink(missing_ok=True)
        _obtain_download(KENT_TOOL_URL, cached, download, reporter, sleep=sleep)
    tool_sha256 = _hash(local)
    reporter.success(
        "kent.validated",
        "Validated the documented bigBedNamedItems command-line signature.",
        pinned_version=KENT_VERSION,
        sha256=tool_sha256,
    )
    return local, tool_sha256


def _archive_hashes(archive: Path, reporter: ProvisioningReporter) -> tuple[str, str]:
    sha256_digest = sha256()
    md5_digest = md5()  # noqa: S324 - UCSC publishes this archive identity as MD5
    total = archive.stat().st_size
    completed = 0
    started_at = time.monotonic()
    with archive.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            sha256_digest.update(chunk)
            md5_digest.update(chunk)
            completed += len(chunk)
            reporter.progress(
                "fasta.archive.verify",
                "Verifying hg38.fa.gz (MD5 + SHA-256)",
                completed,
                total,
                unit="bytes",
                started_at=started_at,
            )
    reporter.progress(
        "fasta.archive.verify",
        "Verifying hg38.fa.gz (MD5 + SHA-256)",
        completed,
        total,
        unit="bytes",
        started_at=started_at,
        force=True,
    )
    return md5_digest.hexdigest(), sha256_digest.hexdigest()


def _gzip_uncompressed_size(archive: Path) -> int:
    with archive.open("rb") as handle:
        handle.seek(-4, os.SEEK_END)
        return int.from_bytes(handle.read(4), "little")


def _decompress_fasta(
    archive: Path,
    destination: Path,
    reporter: ProvisioningReporter,
) -> str:
    temporary = destination.with_name(destination.name + ".part")
    if temporary.exists():
        retained_bytes = temporary.stat().st_size
        reporter.warning(
            "fasta.decompress.restart",
            "A partial decompression was found. Gzip state cannot be resumed safely, so only "
            "this local transform will restart; the verified archive is retained.",
            retained_bytes=retained_bytes,
        )
        temporary.unlink()
    expected_size = _gzip_uncompressed_size(archive)
    free_bytes = shutil.disk_usage(destination.parent).free
    required_bytes = expected_size + 256 * 1024 * 1024
    reporter.info(
        "fasta.disk.preflight",
        f"Drive free space {free_bytes / 1024**3:.2f} GiB; "
        f"decompression requires approximately {required_bytes / 1024**3:.2f} GiB.",
        free_bytes=free_bytes,
        required_bytes=required_bytes,
        expected_fasta_bytes=expected_size,
    )
    if free_bytes < required_bytes:
        raise OSError(
            f"insufficient workspace free space for GRCh38 FASTA: need {required_bytes:,} bytes, "
            f"have {free_bytes:,}"
        )
    digest = sha256()
    completed = 0
    started_at = time.monotonic()
    with gzip.open(archive, "rb") as source, temporary.open("wb") as output:
        while chunk := source.read(4 * 1024 * 1024):
            output.write(chunk)
            digest.update(chunk)
            completed += len(chunk)
            reporter.progress(
                "fasta.decompress",
                "Decompressing hg38.fa.gz",
                completed,
                expected_size,
                unit="bytes",
                started_at=started_at,
            )
        output.flush()
        os.fsync(output.fileno())
    if completed != expected_size:
        raise ValueError(
            f"decompressed FASTA size mismatch: expected {expected_size:,}, observed {completed:,}"
        )
    os.replace(temporary, destination)
    reporter.progress(
        "fasta.decompress",
        "Decompressing hg38.fa.gz",
        completed,
        expected_size,
        unit="bytes",
        started_at=started_at,
        force=True,
    )
    reporter.success(
        "fasta.decompress.complete",
        f"Decompressed hg38.fa ({completed / 1024**3:.2f} GiB).",
        byte_size=completed,
    )
    return digest.hexdigest()


def _prepare_fasta(
    root: Path,
    checkpoint_root: Path,
    download: Download | None,
    reporter: ProvisioningReporter,
    *,
    sleep: Sleeper,
) -> tuple[Path, Path, dict[str, Any]]:
    directory = root / "references/genome/grch38/ucsc-hg38-gca_000001405.15"
    fasta = directory / "hg38.fa"
    index = directory / "hg38.fa.fai"
    completed = directory / "COMPLETED.json"
    ready = directory / "FASTA_READY.json"
    if completed.is_symlink():
        raise ValueError("installed GRCh38 completion manifest is unsafe")
    if completed.is_file():
        try:
            manifest = json.loads(completed.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            manifest = {}
        if (
            manifest.get("schema") == "genome-evidence-reference-completion/v1"
            and manifest.get("source_url") == FASTA_URL
            and manifest.get("upstream_md5") == FASTA_UPSTREAM_MD5
        ):
            if (
                not fasta.is_file()
                or not index.is_file()
                or fasta.is_symlink()
                or index.is_symlink()
            ):
                raise ValueError("completed GRCh38 resource is missing FASTA or FAI")
            reporter.info(
                "fasta.completed.verify",
                "Verifying the existing completed GRCh38 FASTA and FAI.",
                directory=str(directory),
            )
            if _hash(
                fasta,
                reporter=reporter,
                event="fasta.completed.sha256",
            ) != manifest.get("fasta_sha256") or _hash(index) != manifest.get("fai_sha256"):
                raise ValueError("installed GRCh38 FASTA resource failed checksum validation")
            reporter.success(
                "fasta.completed.reuse",
                "Reused the checksum-verified completed GRCh38 FASTA resource.",
                fasta_bytes=fasta.stat().st_size,
            )
            return fasta, index, manifest
        reporter.warning(
            "fasta.completed.repair",
            "The FASTA completion marker is incomplete; verified component checkpoints "
            "will be reused and the marker will be repaired.",
        )
    elif completed.exists():
        raise ValueError("installed GRCh38 completion path is not a regular file")
    directory.mkdir(parents=True, exist_ok=True)
    archive = checkpoint_root / "fasta/ucsc-hg38-gca_000001405.15/hg38.fa.gz"
    partial_archive = archive.with_name(archive.name + ".part")
    if not archive.exists() and partial_archive.is_file():
        reporter.info(
            "fasta.archive.partial",
            f"Found {partial_archive.stat().st_size / 1024**2:.1f} MiB of resumable archive bytes.",
            retained_bytes=partial_archive.stat().st_size,
        )
    _obtain_download(FASTA_URL, archive, download, reporter, sleep=sleep)
    if archive.is_symlink():
        raise ValueError("GRCh38 FASTA archive cache is an unsafe symlink")
    archive_md5, archive_sha256 = _archive_hashes(archive, reporter)
    if archive_md5 != FASTA_UPSTREAM_MD5:
        invalid = archive.with_name(f"{archive.name}.invalid-{archive_sha256[:12]}")
        reporter.warning(
            "fasta.archive.replace",
            "The completed cache file failed UCSC's MD5; preserving it as invalid and "
            "performing one clean download.",
            observed_md5=archive_md5,
            expected_md5=FASTA_UPSTREAM_MD5,
            preserved_path=str(invalid),
        )
        os.replace(archive, invalid)
        _obtain_download(FASTA_URL, archive, download, reporter, sleep=sleep)
        archive_md5, archive_sha256 = _archive_hashes(archive, reporter)
        if archive_md5 != FASTA_UPSTREAM_MD5:
            raise ValueError("GRCh38 FASTA archive does not match UCSC's published MD5")
    fasta_sha256: str | None = None
    if ready.is_symlink():
        raise ValueError("GRCh38 FASTA readiness marker is unsafe")
    if ready.is_file() and fasta.is_file():
        try:
            ready_manifest = json.loads(ready.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            ready_manifest = {}
        if (
            ready_manifest.get("archive_md5") == FASTA_UPSTREAM_MD5
            and ready_manifest.get("byte_size") == fasta.stat().st_size
            and _hash(
                fasta,
                reporter=reporter,
                event="fasta.ready.sha256",
            )
            == ready_manifest.get("fasta_sha256")
        ):
            fasta_sha256 = str(ready_manifest["fasta_sha256"])
            reporter.success(
                "fasta.decompress.resume",
                "Reused the verified decompressed FASTA checkpoint.",
                byte_size=fasta.stat().st_size,
            )
    if fasta_sha256 is None:
        fasta_sha256 = _decompress_fasta(archive, fasta, reporter)
        _json(
            ready,
            {
                "schema": CHECKPOINT_SCHEMA,
                "archive_md5": FASTA_UPSTREAM_MD5,
                "archive_sha256": archive_sha256,
                "fasta_sha256": fasta_sha256,
                "byte_size": fasta.stat().st_size,
            },
        )
    built_index = index.with_name(index.name + ".part")
    build_fasta_index(fasta, built_index, reporter=reporter)
    os.replace(built_index, index)
    reporter.success(
        "fasta.index.complete",
        f"Built and installed hg38.fa.fai ({index.stat().st_size:,} bytes).",
        byte_size=index.stat().st_size,
    )
    manifest = {
        "schema": "genome-evidence-reference-completion/v1",
        "assembly": "GRCh38",
        "assembly_accession": "GCA_000001405.15",
        "source_url": FASTA_URL,
        "upstream_md5": FASTA_UPSTREAM_MD5,
        "archive_sha256": archive_sha256,
        "fasta_sha256": fasta_sha256,
        "fai_sha256": _hash(index),
    }
    _json(completed, manifest)
    try:
        installed_manifest = json.loads(completed.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise OSError("GRCh38 completion marker verification failed") from error
    if installed_manifest != manifest:
        raise OSError("GRCh38 completion marker content mismatch")
    reporter.success(
        "fasta.complete",
        "Published the checksum-verified GRCh38 FASTA, FAI, and completion manifest.",
        directory=str(directory),
    )
    return fasta, index, manifest


def _configured_int(
    environment: Mapping[str, str],
    variable: str,
    default: int,
    minimum: int,
    maximum: int,
) -> int:
    raw = environment.get(variable)
    try:
        value = default if raw is None else int(raw)
    except ValueError as error:
        raise ValueError(f"{variable} must be an integer") from error
    if not minimum <= value <= maximum:
        raise ValueError(f"{variable} must be between {minimum} and {maximum}")
    return value


def _idempotent_json(
    path: Path,
    value: Any,
    reporter: ProvisioningReporter,
    *,
    event: str,
    repair_incomplete: bool,
) -> None:
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if path.is_file():
        if path.is_symlink():
            raise FileExistsError(f"generated artifact path is an unsafe symlink: {path}")
        try:
            installed_payload = path.read_text(encoding="utf-8")
        except UnicodeError as error:
            if not repair_incomplete:
                raise FileExistsError(
                    f"committed generated artifact is not valid UTF-8: {path}"
                ) from error
            installed_payload = None
        if installed_payload == payload:
            reporter.info(f"{event}.resume", f"Reused verified file {path.name}.")
            return
        if not repair_incomplete:
            raise FileExistsError(
                f"existing generated artifact conflicts with expected content: {path}"
            )
        reporter.warning(
            f"{event}.repair",
            f"Repairing an incomplete pre-commit file {path.name}.",
        )
    elif path.exists():
        raise FileExistsError(f"generated artifact path is not a regular file: {path}")
    _json(path, value)
    try:
        installed_payload = path.read_text(encoding="utf-8")
    except UnicodeError as error:
        raise OSError(f"generated JSON verification failed: {path.name}") from error
    if installed_payload != payload:
        raise OSError(f"generated JSON verification failed: {path.name}")
    reporter.success(event, f"Published and verified {path.name}.")


def _idempotent_copy(
    source: Path,
    destination: Path,
    reporter: ProvisioningReporter,
    *,
    event: str,
    repair_incomplete: bool,
) -> None:
    source_sha256 = _hash(source)
    if destination.is_file():
        if destination.is_symlink():
            raise FileExistsError(f"generated artifact path is an unsafe symlink: {destination}")
        if _hash(destination) == source_sha256:
            reporter.info(
                f"{event}.resume",
                f"Reused verified published file {destination.name}.",
                byte_size=destination.stat().st_size,
            )
            return
        if not repair_incomplete:
            raise FileExistsError(
                f"existing generated artifact conflicts with expected content: {destination}"
            )
        reporter.warning(
            f"{event}.repair",
            f"Repairing an incomplete pre-commit file {destination.name}.",
        )
    elif destination.exists():
        raise FileExistsError(f"generated artifact path is not a regular file: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    total = source.stat().st_size
    copied = 0
    started_at = time.monotonic()
    with source.open("rb") as input_handle, temporary.open("wb") as output_handle:
        while chunk := input_handle.read(4 * 1024 * 1024):
            output_handle.write(chunk)
            copied += len(chunk)
            reporter.progress(
                event,
                f"Publishing {destination.name}",
                copied,
                total,
                unit="bytes",
                started_at=started_at,
            )
        output_handle.flush()
        os.fsync(output_handle.fileno())
    if _hash(temporary) != source_sha256:
        raise OSError(f"published copy checksum mismatch: {destination.name}")
    os.replace(temporary, destination)
    if destination.stat().st_size != total or _hash(destination) != source_sha256:
        raise OSError(f"published final checksum mismatch: {destination.name}")
    reporter.progress(
        event,
        f"Publishing {destination.name}",
        copied,
        total,
        unit="bytes",
        started_at=started_at,
        force=True,
    )
    reporter.success(
        f"{event}.complete",
        f"Published and verified {destination.name}.",
        byte_size=destination.stat().st_size,
    )


def _publish_selection(
    root: Path,
    selection: NormalizationResourceSelection,
    *,
    bundle_id: str,
    reporter: ProvisioningReporter,
) -> Path:
    selection_path, publishing_path, completed_path = _selection_control_paths(root)
    value = selection.model_dump(mode="json", by_alias=True)
    expected_payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    staging_path = selection_path.with_name("normalization_resources.staging.json")
    _json(staging_path, value)
    if staging_path.read_text(encoding="utf-8") != expected_payload:
        raise OSError("normalization resource selector staging verification failed")
    selection_sha256 = _hash(staging_path)
    selection_size = staging_path.stat().st_size
    publishing = {
        "schema": SELECTION_PUBLICATION_SCHEMA,
        "source_sha256": selection.source_sha256,
        "bundle_id": bundle_id,
        "selection_sha256": selection_sha256,
        "byte_size": selection_size,
    }
    _json(publishing_path, publishing)
    try:
        installed_publishing = json.loads(publishing_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise OSError("normalization selector publication marker verification failed") from error
    if installed_publishing != publishing:
        raise OSError("normalization selector publication marker content mismatch")
    if selection_path.is_symlink():
        raise FileExistsError("normalization resource selector path is an unsafe symlink")
    if (
        selection_path.is_file()
        and selection_path.stat().st_size == selection_size
        and _hash(selection_path) == selection_sha256
    ):
        reporter.info(
            "selection.publish.resume",
            "Reused the verified selector written by an interrupted publication.",
        )
    else:
        if selection_path.exists() and not selection_path.is_file():
            raise FileExistsError("normalization resource selector path is unsafe")
        _publish_checkpoint_file(staging_path, selection_path)
    completed = {
        "schema": SELECTION_COMPLETION_SCHEMA,
        "source_sha256": selection.source_sha256,
        "bundle_id": bundle_id,
        "sha256": selection_sha256,
        "byte_size": selection_size,
    }
    _json(completed_path, completed)
    if not _selection_completion_matches(root, selection.source_sha256):
        raise OSError("normalization resource selector completion verification failed")
    publishing_path.unlink(missing_ok=True)
    staging_path.unlink(missing_ok=True)
    reporter.success(
        "selection.publish.complete",
        "Published and verified the normalization resource selector and completion marker.",
        selection_path=_workspace_relative(root, selection_path),
    )
    return selection_path


def provision_personal_normalization_resources(
    root: Path,
    environment: Mapping[str, str] | None = None,
    *,
    working_root: Path | None = None,
    download: Download | None = None,
    runner: Runner = subprocess.run,
    sleep: Sleeper = time.sleep,
) -> ProvisioningResult:
    """Provision exact M2 resources with progress, retries, and durable resumption."""
    from .core import validate_workspace

    root = validate_workspace(root)
    env = os.environ if environment is None else environment
    source, source_digest = _selected_source(root, env)
    run_key = sha256(
        f"{source_digest}|dbSNP{DBSNP_BUILD}|Kent{KENT_VERSION}|{FASTA_UPSTREAM_MD5}".encode()
    ).hexdigest()[:24]
    checkpoint_root = root / f"cache/downloads/normalization/v1/{run_key}"
    log_path = root / f"logs/notebooks/00b/{run_key}/events.jsonl"
    reporter = ProvisioningReporter(log_path)
    temporary_parent = working_root.expanduser().resolve() if working_root else None
    if temporary_parent is not None:
        temporary_parent.mkdir(parents=True, exist_ok=True)
    try:
        with reporter:
            checkpoint_root.mkdir(parents=True, exist_ok=True)
            reporter.info(
                "session.start",
                "Starting resumable normalization-resource provisioning. "
                f"Log: {_workspace_relative(root, log_path)} | "
                f"Checkpoint: {_workspace_relative(root, checkpoint_root)}",
                log_path=_workspace_relative(root, log_path),
                checkpoint_path=_workspace_relative(root, checkpoint_root),
                pinned_dbsnp_build=DBSNP_BUILD,
                pinned_kent_version=KENT_VERSION,
            )
            reporter.info(
                "privacy.boundary",
                "Progress contains aggregate counts and public resource metadata only; "
                "genotypes and individual marker identifiers are never logged or transmitted.",
            )
            source_assembly, markers = read_source_markers(
                source, env.get("GENOME_EVIDENCE_SOURCE_BUILD")
            )
            identifiers = tuple(
                sorted(
                    {marker.marker_id for marker in markers if _RSID.fullmatch(marker.marker_id)}
                )
            )
            if not identifiers:
                raise ValueError(
                    "23andMe source contains no canonical numeric rsID markers for dbSNP extraction"
                )
            batch_size = _configured_int(
                env,
                "GENOME_EVIDENCE_DBSNP_BATCH_SIZE",
                _DEFAULT_BATCH_SIZE,
                250,
                25_000,
            )
            query_attempts = _configured_int(
                env,
                "GENOME_EVIDENCE_QUERY_ATTEMPTS",
                _DEFAULT_QUERY_ATTEMPTS,
                1,
                10,
            )
            query_timeout = _configured_int(
                env,
                "GENOME_EVIDENCE_QUERY_TIMEOUT_SECONDS",
                900,
                60,
                3600,
            )
            dbsnp_workers = _configured_int(
                env, "GENOME_EVIDENCE_DBSNP_WORKERS", _DEFAULT_DBSNP_WORKERS, 1, 12
            )
            download_segments = _configured_int(
                env,
                "GENOME_EVIDENCE_COMMON_DOWNLOAD_SEGMENTS",
                _DEFAULT_DOWNLOAD_SEGMENTS,
                1,
                12,
            )
            drive_free = shutil.disk_usage(root).free
            local_free = shutil.disk_usage(temporary_parent or tempfile.gettempdir()).free
            reporter.info(
                "strategy.configured",
                "Strategy: common-first + bounded-parallel fallback. Interrupted work is reusable.",
                configured_workers=dbsnp_workers,
                download_segment_concurrency=download_segments,
                drive_free_bytes=drive_free,
                local_free_bytes=local_free,
            )
            state_path = checkpoint_root / "checkpoint.json"
            if state_path.is_symlink():
                raise ValueError("resource checkpoint identity path is unsafe")
            state: dict[str, Any] | None = None
            if state_path.is_file():
                try:
                    candidate_state = json.loads(state_path.read_text(encoding="utf-8"))
                except (OSError, UnicodeError, json.JSONDecodeError):
                    candidate_state = None
                if (
                    isinstance(candidate_state, dict)
                    and candidate_state.get("schema") == CHECKPOINT_SCHEMA
                ):
                    state = candidate_state
                    if (
                        state.get("source_sha256") != source_digest
                        or state.get("source_assembly") != source_assembly
                        or state.get("dbsnp_build") != DBSNP_BUILD
                        or state.get("kent_version") != KENT_VERSION
                        or state.get("fasta_upstream_md5") != FASTA_UPSTREAM_MD5
                    ):
                        raise ValueError(
                            "resource checkpoint identity conflicts with the selected source"
                        )
                    reporter.info(
                        "checkpoint.resume",
                        "Found a compatible durable checkpoint; completed components will be "
                        "verified and reused.",
                        checkpoint_path=_workspace_relative(root, checkpoint_root),
                    )
                else:
                    reporter.warning(
                        "checkpoint.repair",
                        "The checkpoint identity file is incomplete and will be repaired; "
                        "component completion manifests remain independently verified.",
                    )
            elif state_path.exists():
                raise ValueError("resource checkpoint identity path is not a regular file")
            if state is None:
                state = {
                    "schema": CHECKPOINT_SCHEMA,
                    "source_sha256": source_digest,
                    "source_assembly": source_assembly,
                    "dbsnp_build": DBSNP_BUILD,
                    "kent_version": KENT_VERSION,
                    "fasta_upstream_md5": FASTA_UPSTREAM_MD5,
                    "started_at": datetime.now(UTC).isoformat(),
                }
                _json(state_path, state)
                try:
                    installed_state = json.loads(state_path.read_text(encoding="utf-8"))
                except (OSError, UnicodeError, json.JSONDecodeError) as error:
                    raise OSError("resource checkpoint identity verification failed") from error
                if installed_state != state:
                    raise OSError("resource checkpoint identity content mismatch")
            reporter.success(
                "source.validated",
                f"Validated source metadata: {source_assembly}, {len(markers):,} rows, "
                f"{len(identifiers):,} unique rsIDs.",
                source_assembly=source_assembly,
                source_marker_count=len(markers),
                rsid_count=len(identifiers),
            )
            selection_path, publishing_path, _ = _selection_control_paths(root)
            if publishing_path.is_symlink():
                raise ValueError("normalization selector publication marker is unsafe")
            selection_recovery = publishing_path.exists() and not _selection_completion_matches(
                root, source_digest
            )
            if selection_recovery:
                try:
                    interrupted_publication = json.loads(
                        publishing_path.read_text(encoding="utf-8")
                    )
                except (OSError, UnicodeError, json.JSONDecodeError):
                    interrupted_publication = None
                if (
                    isinstance(interrupted_publication, dict)
                    and interrupted_publication.get("schema") == SELECTION_PUBLICATION_SCHEMA
                    and interrupted_publication.get("source_sha256") != source_digest
                ):
                    raise ValueError(
                        "an interrupted normalization selector publication belongs to a "
                        "different source"
                    )
                reporter.warning(
                    "selection.publish.resume",
                    "Found an interrupted selector publication; verified bundle work will be "
                    "reused and the selector will be repaired at the final step.",
                )
                existing = None
            else:
                existing = load_normalization_resource_selection(
                    root,
                    source_digest,
                    reporter=reporter,
                )
            if existing is not None:
                publishing_path.unlink(missing_ok=True)
                existing_provenance = json.loads(
                    (root / existing.provenance_manifest).read_text(encoding="utf-8")
                )
                statistics = existing_provenance["statistics"]
                defined = int(statistics["defined_marker_count"])
                mapped = (
                    int(statistics["cross_build_mapped_marker_count"])
                    if source_assembly == "GRCh37"
                    else defined
                )
                reporter.success(
                    "session.complete",
                    "All selected resources were already complete and checksum-valid.",
                    defined_marker_count=defined,
                    mapped_marker_count=mapped,
                )
                return ProvisioningResult(
                    existing,
                    selection_path,
                    defined,
                    mapped,
                    max(len(markers) - defined, 0),
                    log_path,
                    checkpoint_root,
                )
            with tempfile.TemporaryDirectory(
                prefix="genome-evidence-resources-", dir=temporary_parent
            ) as tmp:
                work = Path(tmp)
                tool, tool_sha256 = _prepare_kent_tool(
                    root,
                    work,
                    download,
                    runner,
                    reporter,
                    sleep=sleep,
                )
                source_bed, source_query_provenance = _query_common_first(
                    tool,
                    tool_sha256,
                    source_assembly,
                    identifiers,
                    checkpoint_root / "dbsnp",
                    work,
                    runner,
                    reporter,
                    batch_size=batch_size,
                    attempts=query_attempts,
                    timeout_seconds=query_timeout,
                    sleep=sleep,
                    label=source_assembly.lower(),
                    workers=dbsnp_workers,
                    download_segments=download_segments,
                    download=download,
                )
                if source_assembly == "GRCh38":
                    target_bed = source_bed
                    target_query_provenance = source_query_provenance
                else:
                    target_bed, target_query_provenance = _query_common_first(
                        tool,
                        tool_sha256,
                        "GRCh38",
                        identifiers,
                        checkpoint_root / "dbsnp",
                        work,
                        runner,
                        reporter,
                        batch_size=batch_size,
                        attempts=query_attempts,
                        timeout_seconds=query_timeout,
                        sleep=sleep,
                        label="grch38",
                        workers=dbsnp_workers,
                        download_segments=download_segments,
                        download=download,
                    )
                reporter.info(
                    "dbsnp.parse",
                    "Parsing and validating the completed source and target extracts.",
                )
                source_variants = parse_bigbed_variants(source_bed.read_text(encoding="utf-8"))
                target_variants = parse_bigbed_variants(target_bed.read_text(encoding="utf-8"))
                definitions, mappings, stats = build_marker_resources(
                    markers, source_assembly, source_variants, target_variants
                )
                if not definitions:
                    raise ValueError(
                        "dbSNP extraction produced no exact source-build SNV definitions"
                    )
                reporter.success(
                    "markers.built",
                    f"Built {stats['definition_count']:,} exact SNV definitions for "
                    f"{stats['defined_marker_count']:,} source markers.",
                    **stats,
                )
                fasta, fasta_index, fasta_manifest = _prepare_fasta(
                    root,
                    checkpoint_root,
                    download,
                    reporter,
                    sleep=sleep,
                )
                source_extract_sha256 = _hash(source_bed)
                target_extract_sha256 = _hash(target_bed)
                bundle_id = sha256(
                    (
                        f"{source_digest}|dbSNP{DBSNP_BUILD}|{tool_sha256}|"
                        f"{source_extract_sha256}|{target_extract_sha256}|"
                        f"{fasta_manifest['fasta_sha256']}|{RESOURCE_ALGORITHM_VERSION}"
                        f"|{json.dumps(source_query_provenance, sort_keys=True)}"
                        f"|{json.dumps(target_query_provenance, sort_keys=True)}"
                    ).encode()
                ).hexdigest()[:24]
                marker_dir = (
                    root / f"references/markers/23andme/dbsnp{DBSNP_BUILD}-"
                    f"{source_assembly.lower()}-{bundle_id}"
                )
                liftover_dir = (
                    root / f"references/liftover/grch37_to_grch38/dbsnp{DBSNP_BUILD}-{bundle_id}"
                )
                marker_path = marker_dir / "marker-definitions.json"
                source_extract_path = marker_dir / "dbsnp-source-extract.bed"
                target_extract_path = marker_dir / "dbsnp-grch38-extract.bed"
                liftover_path = liftover_dir / "variant-coordinate-map.json"
                provenance_path = root / f"references/manifests/normalization/{bundle_id}.json"
                bundle_completion_path = (
                    root / f"references/manifests/normalization/{bundle_id}.COMPLETED.json"
                )
                bundle_completion = _load_bundle_completion(
                    root,
                    bundle_completion_path,
                    bundle_id=bundle_id,
                    source_sha256=source_digest,
                )
                bundle_committed = bundle_completion is not None
                if bundle_completion_path.exists() and not bundle_committed:
                    reporter.warning(
                        "bundle.publish.repair",
                        "Found an interrupted bundle completion marker; pre-commit files "
                        "will be verified and repaired.",
                    )
                elif bundle_committed:
                    reporter.info(
                        "bundle.publish.resume",
                        "Validated a committed resource bundle; immutable files will be reused.",
                        bundle_id=bundle_id,
                    )
                marker_dir.mkdir(parents=True, exist_ok=True)
                _idempotent_json(
                    marker_path,
                    definitions,
                    reporter,
                    event="publish.marker-definitions",
                    repair_incomplete=not bundle_committed,
                )
                _idempotent_copy(
                    source_bed,
                    source_extract_path,
                    reporter,
                    event="publish.source-extract",
                    repair_incomplete=not bundle_committed,
                )
                _idempotent_copy(
                    target_bed,
                    target_extract_path,
                    reporter,
                    event="publish.target-extract",
                    repair_incomplete=not bundle_committed,
                )
                if source_assembly == "GRCh37":
                    liftover_dir.mkdir(parents=True, exist_ok=True)
                    _idempotent_json(
                        liftover_path,
                        mappings,
                        reporter,
                        event="publish.liftover",
                        repair_incomplete=not bundle_committed,
                    )
                artifacts: dict[str, dict[str, str | int]] = {
                    _workspace_relative(root, marker_path): {
                        "sha256": _hash(marker_path),
                        "byte_size": marker_path.stat().st_size,
                    },
                    _workspace_relative(root, source_extract_path): {
                        "sha256": _hash(source_extract_path),
                        "byte_size": source_extract_path.stat().st_size,
                    },
                    _workspace_relative(root, target_extract_path): {
                        "sha256": _hash(target_extract_path),
                        "byte_size": target_extract_path.stat().st_size,
                    },
                    _workspace_relative(root, fasta): {
                        "sha256": str(fasta_manifest["fasta_sha256"]),
                        "byte_size": fasta.stat().st_size,
                    },
                    _workspace_relative(root, fasta_index): {
                        "sha256": str(fasta_manifest["fai_sha256"]),
                        "byte_size": fasta_index.stat().st_size,
                    },
                }
                if source_assembly == "GRCh37":
                    artifacts[_workspace_relative(root, liftover_path)] = {
                        "sha256": _hash(liftover_path),
                        "byte_size": liftover_path.stat().st_size,
                    }
                provenance: dict[str, Any] = {
                    "schema": PROVENANCE_SCHEMA,
                    "bundle_id": bundle_id,
                    "source_sha256": source_digest,
                    "source_assembly": source_assembly,
                    "retrieved_at": state["started_at"],
                    "vendor_strand_assertion_url": (
                        "https://customercare.23andme.com/hc/en-us/articles/"
                        "212883767-Which-Reference-Genome-and-Strand-Does-23andMe-Use"
                    ),
                    "dbsnp_build": DBSNP_BUILD,
                    "dbsnp_source_url": DBSNP_URLS[source_assembly],
                    "dbsnp_target_url": DBSNP_URLS["GRCh38"],
                    "dbsnp_query_batch_size": batch_size,
                    "dbsnp_workers": dbsnp_workers,
                    "common_download_segments": download_segments,
                    "dbsnp_source_query": source_query_provenance,
                    "dbsnp_target_query": target_query_provenance,
                    "kent_tool_url": KENT_TOOL_URL,
                    "kent_version": KENT_VERSION,
                    "kent_tool_sha256": tool_sha256,
                    "fasta": fasta_manifest,
                    "artifacts": artifacts,
                    "statistics": stats,
                    "checkpoint_schema": CHECKPOINT_SCHEMA,
                    "resource_algorithm_version": RESOURCE_ALGORITHM_VERSION,
                }
                _idempotent_json(
                    provenance_path,
                    provenance,
                    reporter,
                    event="publish.provenance",
                    repair_incomplete=not bundle_committed,
                )
                bundle_completion_value = {
                    "schema": BUNDLE_COMPLETION_SCHEMA,
                    "bundle_id": bundle_id,
                    "source_sha256": source_digest,
                    "algorithm_version": RESOURCE_ALGORITHM_VERSION,
                    "artifacts": artifacts,
                    "provenance": {
                        "path": _workspace_relative(root, provenance_path),
                        "sha256": _hash(provenance_path),
                        "byte_size": provenance_path.stat().st_size,
                    },
                }
                if bundle_completion is not None:
                    if bundle_completion != bundle_completion_value:
                        raise ValueError(
                            "committed normalization resource bundle conflicts with "
                            "expected content"
                        )
                    reporter.info(
                        "bundle.publish.complete.resume",
                        "Reused the verified resource bundle completion marker.",
                    )
                else:
                    _json(bundle_completion_path, bundle_completion_value)
                    try:
                        installed_completion = json.loads(
                            bundle_completion_path.read_text(encoding="utf-8")
                        )
                    except (OSError, UnicodeError, json.JSONDecodeError) as error:
                        raise OSError(
                            "resource bundle completion marker verification failed"
                        ) from error
                    if installed_completion != bundle_completion_value:
                        raise OSError("resource bundle completion marker content mismatch")
                    reporter.success(
                        "bundle.publish.complete",
                        "Published the verified resource bundle completion marker.",
                        bundle_id=bundle_id,
                    )
                selection = NormalizationResourceSelection(
                    source_sha256=source_digest,
                    source_assembly=source_assembly,
                    marker_definitions=_workspace_relative(root, marker_path),
                    grch38_fasta=_workspace_relative(root, fasta),
                    grch37_to_grch38_liftover=(
                        _workspace_relative(root, liftover_path)
                        if source_assembly == "GRCh37"
                        else None
                    ),
                    marker_version=f"UCSC-dbSNP{DBSNP_BUILD}-{source_assembly}",
                    reference_version="UCSC-hg38-GCA_000001405.15",
                    liftover_version=(
                        f"UCSC-dbSNP{DBSNP_BUILD}-matched-rsID"
                        if source_assembly == "GRCh37"
                        else None
                    ),
                    provenance_manifest=_workspace_relative(root, provenance_path),
                )
                selection_path = _publish_selection(
                    root,
                    selection,
                    bundle_id=bundle_id,
                    reporter=reporter,
                )
            defined = int(stats["defined_marker_count"])
            mapped = (
                int(stats["cross_build_mapped_marker_count"])
                if source_assembly == "GRCh37"
                else defined
            )
            reporter.success(
                "session.complete",
                "Provisioning completed; durable selection was published last.",
                defined_marker_count=defined,
                mapped_marker_count=mapped,
                unresolved_marker_count=max(len(markers) - defined, 0),
                selection_path=_workspace_relative(root, selection_path),
            )
            return ProvisioningResult(
                selection,
                selection_path,
                defined,
                mapped,
                max(len(markers) - defined, 0),
                log_path,
                checkpoint_root,
            )
    except FileExistsError:
        raise
    except (_OperationalInterruption, OSError) as error:
        raise ProvisioningIncomplete(
            f"Provisioning is incomplete but resumable: {error}",
            log_path=log_path,
            checkpoint_path=checkpoint_root,
        ) from error
