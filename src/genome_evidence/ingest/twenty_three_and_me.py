"""Source-faithful 23andMe raw-text ingestion (M1 only)."""

import json
import re
import subprocess
from collections import defaultdict
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from typing import Any
from uuid import uuid4

import polars as pl
from pydantic import BaseModel, ConfigDict, Field

from genome_evidence import __version__
from genome_evidence.ingest.base import ParseMode
from genome_evidence.ingest.errors import GenotypeParseError
from genome_evidence.qc.genotype import RECOGNIZED_CHROMOSOMES, categorize_genotype
from genome_evidence.qc.models import (
    AssayQCSummary,
    BuildProvenance,
    CallState,
    ChromosomeQCSummary,
    LexicalGenotypeCategory,
    LexicalZygosity,
    MalformedRecord,
    QCFinding,
    QCSeverity,
    RawGenotypeObservation,
)
from genome_evidence.qc.report import render_qc_report

PARSER_VERSION = "1"
_BUILD_PATTERN = re.compile(r"(?:build|assembly)[\s:=]+(GRCh\d+|hg\d+)", re.IGNORECASE)


class Ingest23andMeConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    mode: ParseMode = ParseMode.STRICT
    genome_build_override: str | None = None
    sample_id: str = Field(default="sample", min_length=1)
    source_identifier: str | None = None
    overwrite: bool = False


class SourceMetadata(BaseModel):
    vendor: str = "23andMe"
    source_format: str = "23andMe raw genotype text"
    logical_source_identifier: str
    source_sha256: str
    file_size_bytes: int
    comment_lines: tuple[str, ...]
    declared_build: str | None
    resolved_build: str
    build_provenance: BuildProvenance
    physical_line_count: int
    parsed_record_count: int
    ingestion_timestamp: datetime
    parser_version: str
    package_version: str
    run_id: str


class IngestionResult(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    run_id: str
    output_directory: Path
    source_metadata: SourceMetadata
    observations: tuple[RawGenotypeObservation, ...]
    malformed_records: tuple[MalformedRecord, ...]
    findings: tuple[QCFinding, ...]
    qc_summary: AssayQCSummary
    manifest: dict[str, Any]


def _git_commit() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, default=str) + "\n").encode()


def _write_bytes(path: Path, content: bytes) -> str:
    path.write_bytes(content)
    return sha256(content).hexdigest()


def _malformed(line_number: int, code: str, explanation: str) -> MalformedRecord:
    return MalformedRecord(
        source_line_number=line_number, reason_code=code, safe_explanation=explanation
    )


def _summarize(
    observations: list[RawGenotypeObservation],
    malformed: list[MalformedRecord],
    findings: list[QCFinding],
    *,
    source_hash: str,
    file_size: int,
    physical: int,
    comments: int,
    blanks: int,
    data: int,
    build: str,
    build_provenance: BuildProvenance,
) -> AssayQCSummary:
    marker_lines: dict[str, list[RawGenotypeObservation]] = defaultdict(list)
    coordinate_lines: dict[tuple[str, int], list[RawGenotypeObservation]] = defaultdict(list)
    exact_lines: dict[tuple[str, str, int, str], list[int]] = defaultdict(list)
    chromosome_records: dict[str, list[RawGenotypeObservation]] = defaultdict(list)
    previous_position: dict[str, int] = {}
    out_of_order = 0
    unknown_chromosomes = 0
    for record in observations:
        marker_lines[record.source_marker_id].append(record)
        coordinate_lines[(record.source_chromosome, record.source_position)].append(record)
        chromosome_records[record.source_chromosome].append(record)
        exact_lines[
            (
                record.source_marker_id,
                record.source_chromosome,
                record.source_position,
                record.raw_genotype,
            )
        ].append(record.source_line_number)
        if record.source_chromosome.upper() not in RECOGNIZED_CHROMOSOMES:
            unknown_chromosomes += 1
            findings.append(
                QCFinding(
                    code="UNRECOGNIZED_CHROMOSOME_TOKEN",
                    severity=QCSeverity.WARNING,
                    message="Source chromosome token is outside the lexical recognition set.",
                    source_line_number=record.source_line_number,
                )
            )
        prior = previous_position.get(record.source_chromosome)
        if prior is not None and record.source_position < prior:
            out_of_order += 1
            findings.append(
                QCFinding(
                    code="OUT_OF_ORDER_RECORD",
                    severity=QCSeverity.INFO,
                    message="Source position decreased within this chromosome token.",
                    source_line_number=record.source_line_number,
                )
            )
        previous_position[record.source_chromosome] = record.source_position
    duplicate_markers = sum(len(group) - 1 for group in marker_lines.values() if len(group) > 1)
    duplicate_coordinates = sum(
        len(group) - 1 for group in coordinate_lines.values() if len(group) > 1
    )
    exact_duplicates = sum(len(lines) - 1 for lines in exact_lines.values() if len(lines) > 1)
    conflicting = 0
    for marker, group in marker_lines.items():
        if len({(r.source_chromosome, r.source_position, r.raw_genotype) for r in group}) > 1:
            conflicting += len(group) - 1
            findings.append(
                QCFinding(
                    code="CONFLICTING_DUPLICATE_MARKER",
                    severity=QCSeverity.WARNING,
                    message=f"Marker {marker!r} has conflicting source records.",
                    related_source_line_numbers=tuple(r.source_line_number for r in group),
                )
            )
        if len(group) > 1:
            findings.append(
                QCFinding(
                    code="DUPLICATE_MARKER_ID",
                    severity=QCSeverity.WARNING,
                    message=f"Marker {marker!r} occurs more than once.",
                    related_source_line_numbers=tuple(r.source_line_number for r in group),
                )
            )
    for group in coordinate_lines.values():
        if len(group) > 1:
            findings.append(
                QCFinding(
                    code="DUPLICATE_COORDINATE",
                    severity=QCSeverity.WARNING,
                    message="A source chromosome/position pair occurs more than once.",
                    related_source_line_numbers=tuple(r.source_line_number for r in group),
                )
            )
    for matching_lines in exact_lines.values():
        if len(matching_lines) > 1:
            findings.append(
                QCFinding(
                    code="EXACT_DUPLICATE_RECORD",
                    severity=QCSeverity.WARNING,
                    message="An exact source record occurs more than once.",
                    related_source_line_numbers=tuple(matching_lines),
                )
            )
    called = sum(r.call_status == CallState.CALLED for r in observations)
    no_calls = len(observations) - called
    chromosome_summaries = []
    for chromosome, group in chromosome_records.items():
        group_called = sum(r.call_status == CallState.CALLED for r in group)
        chromosome_summaries.append(
            ChromosomeQCSummary(
                source_chromosome=chromosome,
                record_count=len(group),
                called_count=group_called,
                no_call_count=len(group) - group_called,
                call_rate=group_called / len(group),
                min_source_position=min(r.source_position for r in group),
                max_source_position=max(r.source_position for r in group),
            )
        )
    conventional = sum(
        r.lexical_category == LexicalGenotypeCategory.TWO_ALLELE_TOKEN
        and r.lexical_zygosity is not None
        for r in observations
    )
    homozygous = sum(r.lexical_zygosity == LexicalZygosity.HOMOZYGOUS_LEXICAL for r in observations)
    heterozygous = sum(
        r.lexical_zygosity == LexicalZygosity.HETEROZYGOUS_LEXICAL for r in observations
    )
    return AssayQCSummary(
        source_sha256=source_hash,
        file_size_bytes=file_size,
        physical_line_count=physical,
        comment_line_count=comments,
        blank_line_count=blanks,
        data_line_count=data,
        parsed_record_count=len(observations),
        malformed_record_count=len(malformed),
        declared_or_resolved_build=build,
        build_provenance=build_provenance,
        eligible_parsed_record_count=len(observations),
        called_record_count=called,
        no_call_record_count=no_calls,
        call_rate=called / len(observations) if observations else None,
        chromosome_summaries=tuple(chromosome_summaries),
        conventional_two_acgt_allele_calls=conventional,
        lexical_homozygous_acgt_calls=homozygous,
        lexical_heterozygous_acgt_calls=heterozygous,
        single_allele_tokens=sum(
            r.lexical_category == LexicalGenotypeCategory.SINGLE_ALLELE_TOKEN for r in observations
        ),
        other_called_tokens=sum(
            r.call_status == CallState.CALLED
            and r.lexical_category == LexicalGenotypeCategory.OTHER_TOKEN
            for r in observations
        ),
        no_calls=no_calls,
        duplicate_marker_id_count=duplicate_markers,
        duplicate_coordinate_count=duplicate_coordinates,
        exact_duplicate_record_count=exact_duplicates,
        conflicting_duplicate_marker_count=conflicting,
        invalid_position_count=sum(r.reason_code == "INVALID_POSITION" for r in malformed),
        unrecognized_chromosome_token_count=unknown_chromosomes,
        out_of_order_record_count=out_of_order,
    )


def _validate_private_output(output_path: Path) -> None:
    """Reject an in-repository output unless Git confirms it is ignored."""
    resolved = output_path.resolve()
    try:
        root = Path(
            subprocess.run(
                ["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True, check=True
            ).stdout.strip()
        ).resolve()
        resolved.relative_to(root)
    except (OSError, subprocess.CalledProcessError, ValueError):
        return
    ignored = (
        subprocess.run(["git", "check-ignore", "-q", str(resolved)], check=False).returncode == 0
    )
    if not ignored:
        raise ValueError("output inside the repository must be covered by .gitignore")


def ingest_23andme(
    input_path: Path, output_path: Path, config: Ingest23andMeConfig | None = None
) -> IngestionResult:
    """Parse one source file and write provenance-linked M1 artifacts.

    Strict mode (the default) aborts at the first malformed data record. Lenient mode
    retains a privacy-safe diagnostic and continues. Neither mode normalizes alleles.
    """
    config = config or Ingest23andMeConfig()
    input_path = input_path.resolve()
    output_path = output_path.resolve()
    if not input_path.is_file():
        raise FileNotFoundError("input source file does not exist")
    _validate_private_output(output_path)
    if output_path.exists() and not output_path.is_dir():
        raise NotADirectoryError("output path exists and is not a directory")
    if output_path.exists() and any(output_path.iterdir()) and not config.overwrite:
        raise FileExistsError("output directory already exists and is not empty")
    output_path.mkdir(parents=True, exist_ok=True)
    raw = input_path.read_bytes()
    source_hash = sha256(raw).hexdigest()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ValueError("source must be UTF-8 text") from error
    lines = text.splitlines()
    run_id = str(uuid4())
    ingested_at = datetime.now(UTC)
    observations: list[RawGenotypeObservation] = []
    malformed: list[MalformedRecord] = []
    comments: list[str] = []
    blank_count = 0
    data_count = 0
    for line_number, line in enumerate(lines, start=1):
        if not line.strip():
            blank_count += 1
            continue
        if line.startswith("#"):
            comments.append(line)
            continue
        data_count += 1
        columns = line.split("\t")
        problem: MalformedRecord | None = None
        position: int | None = None
        if len(columns) != 4:
            problem = _malformed(
                line_number,
                "INVALID_COLUMN_COUNT",
                "Data record must have exactly four tab-separated columns.",
            )
        else:
            try:
                position = int(columns[2])
                if position <= 0:
                    raise ValueError
            except ValueError:
                problem = _malformed(
                    line_number, "INVALID_POSITION", "Source position must be a positive integer."
                )
        if problem is not None:
            if config.mode == ParseMode.STRICT:
                raise GenotypeParseError(problem.source_line_number, problem.safe_explanation)
            malformed.append(problem)
            continue
        assert position is not None
        marker, chromosome, _, genotype = columns
        category, zygosity = categorize_genotype(genotype)
        status = (
            CallState.NO_CALL if category == LexicalGenotypeCategory.NO_CALL else CallState.CALLED
        )
        observations.append(
            RawGenotypeObservation(
                source_marker_id=marker,
                source_chromosome=chromosome,
                source_position=position,
                raw_genotype=genotype,
                call_status=status,
                lexical_category=category,
                lexical_zygosity=zygosity,
                source_line_number=line_number,
                sample_id=config.sample_id,
                ingestion_run_id=run_id,
            )
        )
    declared_build = next(
        (match.group(1) for comment in comments if (match := _BUILD_PATTERN.search(comment))), None
    )
    if config.genome_build_override:
        build, provenance = config.genome_build_override, BuildProvenance.USER_OVERRIDE
    elif declared_build:
        build, provenance = declared_build, BuildProvenance.VENDOR_DECLARED
    else:
        build, provenance = "UNKNOWN", BuildProvenance.UNKNOWN
    findings = [
        QCFinding(
            code="MALFORMED_RECORD",
            severity=QCSeverity.ERROR,
            message=item.safe_explanation,
            source_line_number=item.source_line_number,
        )
        for item in malformed
    ]
    if provenance == BuildProvenance.UNKNOWN:
        findings.append(
            QCFinding(
                code="BUILD_METADATA_ABSENT",
                severity=QCSeverity.INFO,
                message="No explicit genome build metadata was recognized; build remains UNKNOWN.",
            )
        )
    summary = _summarize(
        observations,
        malformed,
        findings,
        source_hash=source_hash,
        file_size=len(raw),
        physical=len(lines),
        comments=len(comments),
        blanks=blank_count,
        data=data_count,
        build=build,
        build_provenance=provenance,
    )
    metadata = SourceMetadata(
        logical_source_identifier=config.source_identifier or f"sha256-{source_hash[:12]}",
        source_sha256=source_hash,
        file_size_bytes=len(raw),
        comment_lines=tuple(comments),
        declared_build=declared_build,
        resolved_build=build,
        build_provenance=provenance,
        physical_line_count=len(lines),
        parsed_record_count=len(observations),
        ingestion_timestamp=ingested_at,
        parser_version=PARSER_VERSION,
        package_version=__version__,
        run_id=run_id,
    )
    config_payload = config.model_dump(mode="json")
    config_hash = sha256(json.dumps(config_payload, sort_keys=True).encode()).hexdigest()
    artifact_hashes: dict[str, str] = {}
    artifact_hashes["source_metadata.json"] = _write_bytes(
        output_path / "source_metadata.json", _json_bytes(metadata.model_dump(mode="json"))
    )
    artifact_hashes["qc_summary.json"] = _write_bytes(
        output_path / "qc_summary.json", _json_bytes(summary.model_dump(mode="json"))
    )
    artifact_hashes["qc_report.md"] = _write_bytes(
        output_path / "qc_report.md", render_qc_report(summary).encode()
    )
    observation_frame = pl.DataFrame(
        [item.model_dump(mode="json") for item in observations],
        schema={
            "source_marker_id": pl.String,
            "source_chromosome": pl.String,
            "source_position": pl.Int64,
            "raw_genotype": pl.String,
            "call_status": pl.String,
            "lexical_category": pl.String,
            "lexical_zygosity": pl.String,
            "source_line_number": pl.Int64,
            "sample_id": pl.String,
            "ingestion_run_id": pl.String,
        },
    )
    observation_frame.write_parquet(output_path / "observations.parquet", compression="zstd")
    artifact_hashes["observations.parquet"] = sha256(
        (output_path / "observations.parquet").read_bytes()
    ).hexdigest()
    finding_frame = pl.DataFrame(
        [item.model_dump(mode="json") for item in findings],
        schema={
            "code": pl.String,
            "severity": pl.String,
            "message": pl.String,
            "source_line_number": pl.Int64,
            "related_source_line_numbers": pl.List(pl.Int64),
        },
    )
    finding_frame.write_parquet(output_path / "qc_findings.parquet", compression="zstd")
    artifact_hashes["qc_findings.parquet"] = sha256(
        (output_path / "qc_findings.parquet").read_bytes()
    ).hexdigest()
    manifest = {
        "run_id": run_id,
        "input_sha256": source_hash,
        "parser_version": PARSER_VERSION,
        "package_version": __version__,
        "git_commit": _git_commit(),
        "configuration": config_payload,
        "configuration_hash": config_hash,
        "started_and_completed_at": ingested_at.isoformat(),
        "artifacts": artifact_hashes,
    }
    _write_bytes(output_path / "manifest.json", _json_bytes(manifest))
    return IngestionResult(
        run_id=run_id,
        output_directory=output_path,
        source_metadata=metadata,
        observations=tuple(observations),
        malformed_records=tuple(malformed),
        findings=tuple(findings),
        qc_summary=summary,
        manifest=manifest,
    )
