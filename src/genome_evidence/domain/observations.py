from pydantic import BaseModel, ConfigDict, Field, model_validator

from genome_evidence.domain.enums import (
    CallStatus,
    LiftoverStatus,
    MappingConfidence,
    MappingStatus,
    StrandTransformation,
)
from genome_evidence.domain.samples import Sample
from genome_evidence.domain.variants import Variant


class GenotypeObservation(BaseModel):
    """Immutable source measurement. Empty alleles encode missingness, never reference."""

    model_config = ConfigDict(frozen=True)
    observation_id: str = Field(min_length=1)
    sample: Sample
    source_marker_id: str = Field(min_length=1)
    source_build: str = Field(min_length=1)
    source_chromosome: str = Field(min_length=1)
    source_position: int = Field(gt=0)
    original_genotype: str | None
    original_alleles: tuple[str, ...] = ()
    call_status: CallStatus
    observation_method: str = Field(min_length=1)
    source_line_identifier: str | None = None

    @model_validator(mode="after")
    def missing_has_no_alleles(self) -> "GenotypeObservation":
        if self.call_status != CallStatus.CALLED and self.original_alleles:
            raise ValueError("missing observations cannot contain called alleles")
        if self.call_status == CallStatus.CALLED and not self.original_alleles:
            raise ValueError("called observations require alleles")
        return self


class ObservationVariantMapping(BaseModel):
    model_config = ConfigDict(frozen=True)
    mapping_id: str = Field(min_length=1)
    observation_id: str = Field(min_length=1)
    variant: Variant | None = None
    status: MappingStatus
    strand_transformation: StrandTransformation
    liftover_status: LiftoverStatus
    is_ambiguous: bool
    confidence: MappingConfidence
    run_id: str = Field(min_length=1)

    @model_validator(mode="after")
    def mapped_requires_variant(self) -> "ObservationVariantMapping":
        if self.status == MappingStatus.MAPPED and self.variant is None:
            raise ValueError("successful mapping requires a canonical variant")
        return self
