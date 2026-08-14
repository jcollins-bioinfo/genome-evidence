"""Typed M2 canonicalization records (derived, never source observations)."""

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, model_validator

from genome_evidence.qc.models import CallState


class MappingOutcome(StrEnum):
    MAPPED = "mapped"
    UNMAPPED = "unmapped"
    AMBIGUOUS = "ambiguous"
    FAILED = "failed"
    UNSUPPORTED = "unsupported"


class StrandTransform(StrEnum):
    NONE = "none"
    COMPLEMENT = "complement"
    REVERSE_COMPLEMENT = "reverse_complement"
    UNKNOWN = "unknown"
    AMBIGUOUS = "ambiguous"


class LiftStatus(StrEnum):
    NOT_REQUIRED = "not_required"
    SUCCESS = "success"
    FAILED = "failed"
    AMBIGUOUS = "ambiguous"
    UNSUPPORTED = "unsupported"


class ReferenceValidation(StrEnum):
    MATCH = "match"
    MISMATCH = "mismatch"
    NOT_ATTEMPTED = "not_attempted"


class ResourceIdentity(BaseModel):
    model_config = ConfigDict(frozen=True)
    resource_type: str
    logical_name: str
    version: str
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    index_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    local_identity: str
    assembly: str | None = None
    source_assembly: str | None = None
    target_assembly: str | None = None
    retrieval_date: str | None = None


class MappingCandidate(BaseModel):
    model_config = ConfigDict(frozen=True)
    candidate_id: str
    observation_reference: str
    normalization_run_id: str
    chromosome: str
    position: int = Field(gt=0)


class CanonicalGenotype(BaseModel):
    model_config = ConfigDict(frozen=True)
    genotype_id: str
    observation_reference: str
    normalization_run_id: str
    variant_id: str
    alleles: tuple[str, ...]
    ploidy: int = Field(ge=1, le=2)
    call_status: CallState

    @model_validator(mode="after")
    def coherent(self) -> "CanonicalGenotype":
        if self.call_status != CallState.CALLED or len(self.alleles) != self.ploidy:
            raise ValueError("canonical genotype must be a called genotype matching ploidy")
        return self


class ObservationMapping(BaseModel):
    model_config = ConfigDict(frozen=True)
    mapping_id: str
    observation_reference: str
    normalization_run_id: str
    outcome: MappingOutcome
    reason: str | None
    source_assembly_token: str
    resolved_source_assembly: str | None
    source_chromosome: str
    source_position: int = Field(gt=0)
    target_assembly: str
    target_chromosome: str | None = None
    target_position: int | None = None
    variant_id: str | None = None
    strand_transform: StrandTransform
    liftover_required: bool
    liftover_status: LiftStatus
    reference_validation: ReferenceValidation
    candidate_ids: tuple[str, ...] = ()
    source_to_canonical_alleles: tuple[str, ...] = ()

    @model_validator(mode="after")
    def coherent(self) -> "ObservationMapping":
        if self.outcome == MappingOutcome.MAPPED:
            if not all((self.variant_id, self.target_chromosome, self.target_position)):
                raise ValueError("mapped outcome requires a canonical target")
            if self.reference_validation != ReferenceValidation.MATCH:
                raise ValueError("mapped outcome requires successful REF validation")
            if self.liftover_required and self.liftover_status != LiftStatus.SUCCESS:
                raise ValueError("mapped lifted outcome requires successful liftover")
        elif self.variant_id is not None:
            raise ValueError("non-mapped outcome cannot claim a canonical variant")
        if self.outcome == MappingOutcome.AMBIGUOUS and not self.candidate_ids:
            raise ValueError("ambiguous outcome requires candidates")
        return self
