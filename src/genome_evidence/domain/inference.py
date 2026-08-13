from pydantic import BaseModel, ConfigDict, Field, model_validator

from genome_evidence.domain.enums import PhasingStatus
from genome_evidence.domain.variants import Variant


class GenotypeInference(BaseModel):
    """A statistical genotype result, never a direct observation."""

    model_config = ConfigDict(frozen=True)
    inference_id: str = Field(min_length=1)
    sample_id: str = Field(min_length=1)
    variant: Variant
    genotype_probabilities: dict[str, float] = Field(default_factory=dict)
    dosage: float | None = Field(default=None, ge=0, le=2)
    imputation_quality: float | None = Field(default=None, ge=0, le=1)
    reference_panel: str = Field(min_length=1)
    reference_panel_version: str = Field(min_length=1)
    phasing_status: PhasingStatus
    run_id: str = Field(min_length=1)

    @model_validator(mode="after")
    def probabilities_are_valid(self) -> "GenotypeInference":
        if any(value < 0 or value > 1 for value in self.genotype_probabilities.values()):
            raise ValueError("genotype probabilities must be between zero and one")
        total = sum(self.genotype_probabilities.values())
        if self.genotype_probabilities and abs(total - 1.0) > 1e-6:
            raise ValueError("genotype probabilities must sum to one")
        return self
