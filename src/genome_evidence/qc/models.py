"""Typed descriptive assay-QC results; none are biological interpretations."""

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field


class BuildProvenance(StrEnum):
    VENDOR_DECLARED = "vendor_declared"
    USER_OVERRIDE = "user_override"
    UNKNOWN = "unknown"


class CallState(StrEnum):
    CALLED = "called"
    NO_CALL = "no_call"
    MALFORMED = "malformed"


class LexicalGenotypeCategory(StrEnum):
    NO_CALL = "no_call"
    SINGLE_ALLELE_TOKEN = "single_allele_token"
    TWO_ALLELE_TOKEN = "two_allele_token"
    OTHER_TOKEN = "other_token"


class LexicalZygosity(StrEnum):
    HOMOZYGOUS_LEXICAL = "homozygous_lexical"
    HETEROZYGOUS_LEXICAL = "heterozygous_lexical"


class QCSeverity(StrEnum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


class RawGenotypeObservation(BaseModel):
    model_config = ConfigDict(frozen=True)
    source_marker_id: str
    source_chromosome: str
    source_position: int = Field(gt=0)
    raw_genotype: str
    call_status: CallState
    lexical_category: LexicalGenotypeCategory
    lexical_zygosity: LexicalZygosity | None = None
    source_line_number: int = Field(gt=0)
    sample_id: str
    ingestion_run_id: str


class MalformedRecord(BaseModel):
    model_config = ConfigDict(frozen=True)
    source_line_number: int = Field(gt=0)
    reason_code: str
    safe_explanation: str


class QCFinding(BaseModel):
    model_config = ConfigDict(frozen=True)
    code: str
    severity: QCSeverity
    message: str
    source_line_number: int | None = None
    related_source_line_numbers: tuple[int, ...] = ()


class ChromosomeQCSummary(BaseModel):
    source_chromosome: str
    record_count: int
    called_count: int
    no_call_count: int
    call_rate: float | None
    min_source_position: int
    max_source_position: int


class AssayQCSummary(BaseModel):
    source_sha256: str
    file_size_bytes: int
    physical_line_count: int
    comment_line_count: int
    blank_line_count: int
    data_line_count: int
    parsed_record_count: int
    malformed_record_count: int
    declared_or_resolved_build: str
    build_provenance: BuildProvenance
    eligible_parsed_record_count: int
    called_record_count: int
    no_call_record_count: int
    call_rate: float | None
    chromosome_summaries: tuple[ChromosomeQCSummary, ...]
    conventional_two_acgt_allele_calls: int
    lexical_homozygous_acgt_calls: int
    lexical_heterozygous_acgt_calls: int
    single_allele_tokens: int
    other_called_tokens: int
    no_calls: int
    duplicate_marker_id_count: int
    duplicate_coordinate_count: int
    exact_duplicate_record_count: int
    conflicting_duplicate_marker_count: int
    invalid_position_count: int
    unrecognized_chromosome_token_count: int
    out_of_order_record_count: int
