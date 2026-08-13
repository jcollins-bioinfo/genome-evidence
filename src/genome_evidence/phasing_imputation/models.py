"""Immutable M6 contracts; all outputs are explicitly model inferences."""

from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, model_validator


class FrozenModel(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


class Artifact(FrozenModel):
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    byte_size: int = Field(ge=0)


class ChromosomeReference(FrozenModel):
    chromosome: str
    panel: str
    panel_index: str
    genetic_map: str
    variants: str
    artifacts: dict[str, Artifact]


class ReferenceManifest(FrozenModel):
    schema_id: str = Field(alias="schema")
    bundle_id: str
    assembly: str
    chromosomes: tuple[str, ...]
    source_release: str
    source_urls: tuple[str, ...]
    citations: tuple[str, ...]
    license: str
    redistribution_restrictions: str
    phased_diploid: bool
    normalization_fasta_identity: str
    map_provenance: str
    map_units: str
    allele_frequency_semantics: str
    sample_qc: str
    variant_qc: str
    commands: tuple[str, ...]
    software: tuple[str, ...]
    intended_use: str
    limitations: str
    engine_compatibility: tuple[str, ...]
    per_chromosome: tuple[ChromosomeReference, ...]

    @model_validator(mode="after")
    def supported(self) -> "ReferenceManifest":
        expected = tuple(str(i) for i in range(1, 23))
        if self.schema_id != "genome-evidence-phasing-reference/v1" or self.assembly != "GRCh38":
            raise ValueError("unsupported phasing reference schema or assembly")
        if not self.phased_diploid or self.map_units != "cM":
            raise ValueError("reference must be phased diploid with centimorgan maps")
        if not self.chromosomes or any(x not in expected for x in self.chromosomes):
            raise ValueError("reference scope must contain only autosomes 1-22")
        if tuple(x.chromosome for x in self.per_chromosome) != self.chromosomes:
            raise ValueError("per-chromosome reference declarations are incomplete or unordered")
        return self


class M6Config(FrozenModel):
    chromosomes: tuple[str, ...] = ("22",)
    seed: int = 20260101
    threads: int = Field(default=1, ge=1)
    memory_mb: int = Field(default=2048, ge=512)
    ne: int = Field(default=1_000_000, gt=0)
    timeout_seconds: int = Field(default=3600, gt=0)
    masked_fraction: float = Field(default=0.05, ge=0, lt=1)

    @model_validator(mode="after")
    def autosomes(self) -> "M6Config":
        if not self.chromosomes or len(set(self.chromosomes)) != len(self.chromosomes):
            raise ValueError("chromosome scope must be unique and nonempty")
        if any(c not in {str(i) for i in range(1, 23)} for c in self.chromosomes):
            raise ValueError("M6 supports chromosomes 1-22 only")
        return self


class RunState(StrEnum):
    PLANNED = "planned"
    ALIGNED = "aligned"
    TARGET_WRITTEN = "target_written"
    ENGINE_COMPLETED = "engine_completed"
    PARSED = "parsed"
    VALIDATED = "validated"


class ReferenceValidation(FrozenModel):
    directory: Path
    manifest: ReferenceManifest
    manifest_sha256: str
