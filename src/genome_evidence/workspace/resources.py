"""Pinned, provenance-bearing normalization-resource provisioning for personal data."""

import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from hashlib import md5, sha256
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from genome_evidence.normalization.resources import canonical_assembly, canonical_chromosome

DBSNP_BUILD = "155"
DBSNP_URLS = {
    "GRCh37": "https://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp155.bb",
    "GRCh38": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155.bb",
}
FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
FASTA_UPSTREAM_MD5 = "1c9dcaddfa41027f17cd8f7a82c7293b"
KENT_VERSION = "479"
KENT_TOOL_URL = (
    f"https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v{KENT_VERSION}/bigBedNamedItems"
)
SELECTION_SCHEMA = "genome-evidence-normalization-resource-selection/v1"
PROVENANCE_SCHEMA = "genome-evidence-normalization-resource-provenance/v1"
_BUILD_PATTERN = re.compile(r"(?:build|assembly)[\s:=]+(GRCh\d+|hg\d+|37|38)\b", re.IGNORECASE)
_DIGEST = re.compile(r"^[0-9a-f]{64}$")


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


def load_normalization_resource_selection(
    root: Path, source_sha256: str
) -> NormalizationResourceSelection | None:
    """Load a durable selection and prove that it belongs to the selected private source."""
    path = root / "config/normalization_resources.json"
    if not path.exists():
        return None
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
    artifact_paths = [
        selection.marker_definitions,
        selection.grch38_fasta,
    ]
    if selection.grch37_to_grch38_liftover is not None:
        artifact_paths.append(selection.grch37_to_grch38_liftover)
    relative_paths = [*artifact_paths, selection.provenance_manifest]
    resolved: dict[str, Path] = {}
    for relative in relative_paths:
        candidate = Path(relative)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError("normalization resource selection paths must be workspace-relative")
        resource = (root / candidate).resolve()
        if (
            not resource.is_relative_to(root.resolve())
            or not resource.is_file()
            or resource.is_symlink()
        ):
            raise ValueError("normalization resource selection references a missing or unsafe file")
        resolved[relative] = resource
    try:
        provenance = json.loads(resolved[selection.provenance_manifest].read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
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
    for relative in artifact_paths:
        identity = artifacts.get(relative)
        if not isinstance(identity, dict) or not _DIGEST.fullmatch(str(identity.get("sha256"))):
            raise ValueError("normalization resource provenance does not cover every selection")
        resource = resolved[relative]
        if resource.stat().st_size != identity.get("byte_size"):
            raise ValueError("normalization resource size does not match provenance")
        if _hash(resource) != identity["sha256"]:
            raise ValueError("normalization resource checksum does not match provenance")
    return selection


Download = Callable[[str, Path], None]
Runner = Callable[..., subprocess.CompletedProcess[str]]


def _hash(path: Path, algorithm: str = "sha256") -> str:
    digest = sha256() if algorithm == "sha256" else md5()  # noqa: S324 - upstream publishes MD5
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "genome-evidence/0.4"})
    with urllib.request.urlopen(request, timeout=120) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)


def _json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _workspace_relative(root: Path, path: Path) -> str:
    resolved = path.resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError("resource path escapes the private workspace")
    return str(resolved.relative_to(root.resolve()))


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
        except (OSError, json.JSONDecodeError):
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


def build_fasta_index(fasta: Path, output: Path) -> None:
    """Build a standard .fai while rejecting irregular non-terminal sequence lines."""
    entries: list[tuple[str, int, int, int, int]] = []
    name: str | None = None
    length = 0
    offset = 0
    line_bases = 0
    line_width = 0
    previous_bases: int | None = None
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


def _run_bigbed_query(
    tool: Path,
    url: str,
    identifiers: Path,
    output: Path,
    runner: Runner,
) -> None:
    runner(
        [str(tool), "-nameFile", url, str(identifiers), str(output)],
        check=True,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if not output.is_file():
        raise RuntimeError("UCSC query completed without producing an output file")


def _prepare_kent_tool(
    root: Path, work: Path, download: Download, runner: Runner
) -> tuple[Path, str]:
    cached = root / f"cache/tools/ucsc/kent-v{KENT_VERSION}/bigBedNamedItems"
    if not cached.is_file():
        temporary = work / "bigBedNamedItems.download"
        download(KENT_TOOL_URL, temporary)
        os.chmod(temporary, 0o700)
        cached.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(temporary, cached)
        os.chmod(cached, 0o700)
    local = work / "bigBedNamedItems"
    shutil.copyfile(cached, local)
    os.chmod(local, 0o700)
    probe = runner([str(local)], capture_output=True, text=True, check=False, timeout=30)
    probe_text = (probe.stdout or "") + (probe.stderr or "")
    if f"kent source version {KENT_VERSION}" not in probe_text.lower():
        raise ValueError("cached UCSC utility does not report the pinned Kent version")
    return local, _hash(local)


def _prepare_fasta(root: Path, work: Path, download: Download) -> tuple[Path, Path, dict[str, Any]]:
    directory = root / "references/genome/grch38/ucsc-hg38-gca_000001405.15"
    fasta = directory / "hg38.fa"
    index = directory / "hg38.fa.fai"
    completed = directory / "COMPLETED.json"
    if completed.is_file() and fasta.is_file() and index.is_file():
        manifest = json.loads(completed.read_text(encoding="utf-8"))
        if _hash(fasta) != manifest.get("fasta_sha256") or _hash(index) != manifest.get(
            "fai_sha256"
        ):
            raise ValueError("installed GRCh38 FASTA resource failed checksum validation")
        return fasta, index, manifest
    if directory.exists():
        raise FileExistsError("incomplete GRCh38 resource directory exists; preserve or remove it")
    archive = work / "hg38.fa.gz"
    unpacked = work / "hg38.fa"
    built_index = work / "hg38.fa.fai"
    download(FASTA_URL, archive)
    if _hash(archive, "md5") != FASTA_UPSTREAM_MD5:
        raise ValueError("GRCh38 FASTA archive does not match UCSC's published MD5")
    with gzip.open(archive, "rb") as source, unpacked.open("wb") as destination:
        shutil.copyfileobj(source, destination, length=1024 * 1024)
    build_fasta_index(unpacked, built_index)
    directory.mkdir(parents=True)
    shutil.copyfile(unpacked, fasta)
    shutil.copyfile(built_index, index)
    manifest = {
        "schema": "genome-evidence-reference-completion/v1",
        "assembly": "GRCh38",
        "assembly_accession": "GCA_000001405.15",
        "source_url": FASTA_URL,
        "upstream_md5": FASTA_UPSTREAM_MD5,
        "archive_sha256": _hash(archive),
        "fasta_sha256": _hash(fasta),
        "fai_sha256": _hash(index),
    }
    _json(completed, manifest)
    return fasta, index, manifest


def provision_personal_normalization_resources(
    root: Path,
    environment: Mapping[str, str] | None = None,
    *,
    working_root: Path | None = None,
    download: Download = _download,
    runner: Runner = subprocess.run,
) -> ProvisioningResult:
    """Provision the exact source-compatible M2 resources into canonical Drive subdirectories."""
    from .core import validate_workspace

    root = validate_workspace(root)
    env = os.environ if environment is None else environment
    source, source_digest = _selected_source(root, env)
    source_assembly, markers = read_source_markers(source, env.get("GENOME_EVIDENCE_SOURCE_BUILD"))
    existing = load_normalization_resource_selection(root, source_digest)
    if existing is not None:
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
        return ProvisioningResult(
            existing,
            root / "config/normalization_resources.json",
            defined,
            mapped,
            len(markers) - defined,
        )
    temporary_parent = working_root.expanduser().resolve() if working_root else None
    if temporary_parent is not None:
        temporary_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="genome-evidence-resources-", dir=temporary_parent
    ) as tmp:
        work = Path(tmp)
        identifiers = work / "rsids.txt"
        identifiers.write_text(
            "\n".join(sorted({x.marker_id for x in markers if x.marker_id.startswith("rs")}))
            + "\n",
            encoding="ascii",
        )
        tool, tool_sha256 = _prepare_kent_tool(root, work, download, runner)
        source_bed = work / "dbsnp-source.bed"
        target_bed = work / "dbsnp-target.bed"
        _run_bigbed_query(tool, DBSNP_URLS[source_assembly], identifiers, source_bed, runner)
        if source_assembly == "GRCh38":
            shutil.copyfile(source_bed, target_bed)
        else:
            _run_bigbed_query(tool, DBSNP_URLS["GRCh38"], identifiers, target_bed, runner)
        source_variants = parse_bigbed_variants(source_bed.read_text(encoding="utf-8"))
        target_variants = parse_bigbed_variants(target_bed.read_text(encoding="utf-8"))
        definitions, mappings, stats = build_marker_resources(
            markers, source_assembly, source_variants, target_variants
        )
        if not definitions:
            raise ValueError("dbSNP extraction produced no exact source-build SNV definitions")
        fasta, fasta_index, fasta_manifest = _prepare_fasta(root, work, download)
        bundle_id = sha256(
            f"{source_digest}|dbSNP{DBSNP_BUILD}|{fasta_manifest['fasta_sha256']}".encode()
        ).hexdigest()[:24]
        marker_dir = (
            root
            / f"references/markers/23andme/dbsnp{DBSNP_BUILD}-{source_assembly.lower()}-{bundle_id}"
        )
        liftover_dir = root / f"references/liftover/grch37_to_grch38/dbsnp{DBSNP_BUILD}-{bundle_id}"
        marker_path = marker_dir / "marker-definitions.json"
        source_extract_path = marker_dir / "dbsnp-source-extract.bed"
        target_extract_path = marker_dir / "dbsnp-grch38-extract.bed"
        liftover_path = liftover_dir / "variant-coordinate-map.json"
        if marker_dir.exists():
            raise FileExistsError("marker resource destination already exists")
        marker_dir.mkdir(parents=True)
        _json(marker_path, definitions)
        shutil.copyfile(source_bed, source_extract_path)
        shutil.copyfile(target_bed, target_extract_path)
        if source_assembly == "GRCh37":
            liftover_dir.mkdir(parents=True)
            _json(liftover_path, mappings)
        retrieved = datetime.now(UTC).isoformat()
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
                "sha256": _hash(fasta),
                "byte_size": fasta.stat().st_size,
            },
            _workspace_relative(root, fasta_index): {
                "sha256": _hash(fasta_index),
                "byte_size": fasta_index.stat().st_size,
            },
        }
        provenance: dict[str, Any] = {
            "schema": PROVENANCE_SCHEMA,
            "bundle_id": bundle_id,
            "source_sha256": source_digest,
            "source_assembly": source_assembly,
            "retrieved_at": retrieved,
            "vendor_strand_assertion_url": (
                "https://customercare.23andme.com/hc/en-us/articles/"
                "212883767-Which-Reference-Genome-and-Strand-Does-23andMe-Use"
            ),
            "dbsnp_build": DBSNP_BUILD,
            "dbsnp_source_url": DBSNP_URLS[source_assembly],
            "dbsnp_target_url": DBSNP_URLS["GRCh38"],
            "kent_tool_url": KENT_TOOL_URL,
            "kent_version": KENT_VERSION,
            "kent_tool_sha256": tool_sha256,
            "fasta": fasta_manifest,
            "artifacts": artifacts,
            "statistics": stats,
        }
        if source_assembly == "GRCh37":
            artifacts[_workspace_relative(root, liftover_path)] = {
                "sha256": _hash(liftover_path),
                "byte_size": liftover_path.stat().st_size,
            }
        provenance_path = root / f"references/manifests/normalization/{bundle_id}.json"
        _json(provenance_path, provenance)
        selection = NormalizationResourceSelection(
            source_sha256=source_digest,
            source_assembly=source_assembly,
            marker_definitions=_workspace_relative(root, marker_path),
            grch38_fasta=_workspace_relative(root, fasta),
            grch37_to_grch38_liftover=(
                _workspace_relative(root, liftover_path) if source_assembly == "GRCh37" else None
            ),
            marker_version=f"UCSC-dbSNP{DBSNP_BUILD}-{source_assembly}",
            reference_version="UCSC-hg38-GCA_000001405.15",
            liftover_version=(
                f"UCSC-dbSNP{DBSNP_BUILD}-matched-rsID" if source_assembly == "GRCh37" else None
            ),
            provenance_manifest=_workspace_relative(root, provenance_path),
        )
        selection_path = root / "config/normalization_resources.json"
        _json(selection_path, selection.model_dump(mode="json", by_alias=True))
    defined = stats["defined_marker_count"]
    mapped = stats["cross_build_mapped_marker_count"] if source_assembly == "GRCh37" else defined
    return ProvisioningResult(selection, selection_path, defined, mapped, len(markers) - defined)
