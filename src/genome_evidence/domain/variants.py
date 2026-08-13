from pydantic import BaseModel, ConfigDict, Field, field_validator


class Variant(BaseModel):
    """A canonical allele; rsID is metadata, not identity."""

    model_config = ConfigDict(frozen=True)
    assembly: str = Field(min_length=1)
    chromosome: str = Field(min_length=1)
    position: int = Field(gt=0)
    reference: str = Field(min_length=1)
    alternate: str = Field(min_length=1)
    rsid: str | None = None

    @field_validator("reference", "alternate")
    @classmethod
    def normalize_allele_case(cls, value: str) -> str:
        return value.upper()
