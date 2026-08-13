"""Typed, source-faithful M3 external assertion and linking boundaries."""

from datetime import date, datetime
from enum import StrEnum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator


class AssertionLevel(StrEnum):
    SUBMITTED = "source_submitted"
    AGGREGATE = "source_computed_aggregate"


class ClassificationType(StrEnum):
    GERMLINE = "germline"
    SOMATIC_CLINICAL_IMPACT = "somatic_clinical_impact"
    ONCOGENICITY = "oncogenicity"
    UNKNOWN = "unknown"


class SourceRecordStatus(StrEnum):
    CURRENT = "current"
    REPLACED = "replaced"
    DELETED = "deleted"
    UNKNOWN = "unknown"


class LinkOutcome(StrEnum):
    MATCHED = "matched"
    UNMATCHED = "unmatched"
    AMBIGUOUS = "ambiguous"
    INCOMPATIBLE = "incompatible"
    UNSUPPORTED = "unsupported"


class LinkReason(StrEnum):
    EXACT_ALLELE_MATCH = "exact_allele_match"
    ALLELE_NOT_IN_NORMALIZATION_RUN = "allele_not_in_normalization_run"
    MULTIPLE_CANONICAL_TARGETS = "multiple_canonical_targets"
    ASSEMBLY_INCOMPATIBLE = "assembly_incompatible"
    UNSUPPORTED_VARIATION_TYPE = "unsupported_variation_type"
    INCOMPLETE_ALLELE_IDENTITY = "incomplete_allele_identity"


class ExternalSourceSnapshot(BaseModel):
    model_config = ConfigDict(frozen=True)
    snapshot_id: str
    source_namespace: str
    dataset: str
    release_identity: str
    release_date: date
    retrieved_at: datetime
    source_file_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    source_file_size_bytes: int = Field(ge=0)
    xml_format: str
    xml_schema_identity: str | None = None
    release_identity_override: str | None = None


class ExternalVariantRepresentation(BaseModel):
    model_config = ConfigDict(frozen=True)
    representation_id: str
    snapshot_id: str
    vcv_accession: str
    vcv_version: int = Field(ge=1)
    variation_id: int | None = None
    allele_id: int | None = None
    variation_type: str
    assembly: str | None = None
    chromosome: str | None = None
    position: int | None = Field(default=None, gt=0)
    reference: str | None = None
    alternate: str | None = None
    source_rsid: str | None = None


class ExternalAssertion(BaseModel):
    model_config = ConfigDict(frozen=True)
    assertion_instance_id: str
    logical_source_key: str
    snapshot_id: str
    representation_id: str
    source_accession: str
    accession_version: int = Field(ge=1)
    vcv_accession: str
    rcv_accessions: tuple[str, ...] = ()
    scv_accession: str | None = None
    assertion_level: AssertionLevel
    classification_type: ClassificationType
    source_classification_terms: tuple[str, ...]
    source_review_status: str | None = None
    source_record_status: SourceRecordStatus = SourceRecordStatus.UNKNOWN
    source_record_status_term: str | None = None
    submitter_name: str | None = None
    submitter_identifier: str | None = None
    date_last_evaluated: date | None = None
    assertion_method: str | None = None
    citation_identifiers: tuple[str, ...] = ()
    observed_in_count: int = 0
    source_evidence_structure_count: int = 0
    condition_ids: tuple[str, ...] = ()
    record_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    extension_payload: dict[str, Any] = Field(default_factory=dict)


class AssertionRelationship(BaseModel):
    model_config = ConfigDict(frozen=True)
    relationship_id: str
    snapshot_id: str
    subject_assertion_id: str
    predicate: str
    object_assertion_id: str


class ExternalCondition(BaseModel):
    model_config = ConfigDict(frozen=True)
    condition_id: str
    snapshot_id: str
    source_identifier: str | None = None
    source_identifier_type: str | None = None
    source_name: str


class AssertionConditionRelationship(BaseModel):
    model_config = ConfigDict(frozen=True)
    relationship_id: str
    assertion_instance_id: str
    condition_id: str


class VariantEvidenceLink(BaseModel):
    model_config = ConfigDict(frozen=True)
    link_id: str
    annotation_run_id: str
    representation_id: str
    variant_id: str | None = None
    outcome: LinkOutcome
    reason: LinkReason
    assembly: str | None = None
    chromosome: str | None = None
    position: int | None = None
    reference: str | None = None
    alternate: str | None = None

    @model_validator(mode="after")
    def matched_is_exact(self) -> "VariantEvidenceLink":
        if self.outcome == LinkOutcome.MATCHED and not self.variant_id:
            raise ValueError("a matched evidence link requires an M2 variant ID")
        if self.outcome != LinkOutcome.MATCHED and self.variant_id is not None:
            raise ValueError("only a matched evidence link may identify an M2 variant")
        return self


class ClinVarIngestionConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    release_identity_override: str | None = None
    retrieval_timestamp: datetime | None = None


class EvidenceIngestionResult(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    run_id: str
    output_directory: Path
    snapshot: ExternalSourceSnapshot
    variants: tuple[ExternalVariantRepresentation, ...]
    assertions: tuple[ExternalAssertion, ...]
    relationships: tuple[AssertionRelationship, ...]
    conditions: tuple[ExternalCondition, ...]
    manifest: dict[str, Any]


class AnnotationConfig(BaseModel):
    model_config = ConfigDict(frozen=True)


class AnnotationResult(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    run_id: str
    output_directory: Path
    links: tuple[VariantEvidenceLink, ...]
    manifest: dict[str, Any]
