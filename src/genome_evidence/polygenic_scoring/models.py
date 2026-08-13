"""Typed, immutable contracts for narrow additive polygenic scoring."""

from decimal import Decimal
from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, field_validator


class SourcePolicy(StrEnum):
    OBSERVED_ONLY = "observed_only"
    OBSERVED_THEN_IMPUTED = "observed_then_imputed"


class Evaluability(StrEnum):
    EVALUABLE = "evaluable"
    PARTIAL_NOT_EVALUABLE = "partial_not_evaluable"
    NOT_EVALUABLE = "not_evaluable"


class ScoreConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    pgs_ids: tuple[str, ...] = Field(min_length=1)
    genotype_source_policy: SourcePolicy = SourcePolicy.OBSERVED_ONLY
    reference_distribution: str | None = None

    @field_validator("pgs_ids")
    @classmethod
    def explicit_unique_ids(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        if len(set(value)) != len(value) or any(not x.startswith("PGS") for x in value):
            raise ValueError("PGS IDs must be explicit, unique PGS Catalog identifiers")
        return value


class BundleValidation(BaseModel):
    model_config = ConfigDict(frozen=True)
    directory: Path
    bundle_id: str
    model_ids: tuple[str, ...]
    artifact_count: int


class ScoreRunResult(BaseModel):
    model_config = ConfigDict(frozen=True)
    run_id: str
    output_directory: Path
    model_count: int
    statuses: dict[str, Evaluability]


def decimal_string(value: Decimal) -> str:
    """Serialize Decimal deterministically without binary floating-point conversion."""
    return format(value.normalize(), "f") if value else "0"
