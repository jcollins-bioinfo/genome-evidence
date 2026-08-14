"""Strict M8 pharmacogenomic reference and candidate-evidence types.

These frozen models keep nomenclature, source assertions, observed-locus evidence,
software candidates, and phenotype evidence separate.  They deliberately cannot
represent a prescription or clinical laboratory result.
"""

from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field


class GeneStrategy(StrEnum):
    """Declared evaluation capability for a pharmacogene."""

    STAR_HAPLOTYPE_SMALL_VARIANT = "star_haplotype_small_variant"
    SINGLE_VARIANT_EVIDENCE = "single_variant_evidence"
    EXTERNAL_VALIDATED_CALL_REQUIRED = "external_validated_call_required"
    UNSUPPORTED = "unsupported"


class ConstraintState(StrEnum):
    """Explicit meaning of one source-normalized allele-definition cell."""

    REQUIRED_REFERENCE = "required_reference"
    REQUIRED_ALTERNATE = "required_alternate"
    ALLOWED_ALLELES = "allowed_alleles"
    NOT_CONSTRAINING = "not_constraining"
    UNSUPPORTED_STRUCTURAL_REQUIREMENT = "unsupported_structural_requirement"


class LocusEvidenceState(StrEnum):
    """Relationship between an exact bundle locus and M2 observations."""

    OBSERVED_REFERENCE = "observed_reference"
    OBSERVED_HETEROZYGOUS = "observed_heterozygous"
    OBSERVED_ALTERNATE = "observed_alternate"
    MISSING_OR_UNASSAYED = "missing_or_unassayed"
    NO_CALL = "no_call"
    DUPLICATE_CONCORDANT = "duplicate_concordant"
    DUPLICATE_CONFLICT = "duplicate_conflict"
    ALLELE_INCOMPATIBLE = "allele_incompatible"
    PLOIDY_UNSUPPORTED = "ploidy_unsupported"
    MAPPING_UNRESOLVED = "mapping_unresolved"
    UNSUPPORTED_VARIANT_CLASS = "unsupported_variant_class"
    UNMODELED_OBSERVED_VARIANT = "unmodeled_observed_variant"


class GeneOutcome(StrEnum):
    """Fail-closed outcome of generic candidate enumeration."""

    RESOLVED_CANDIDATE = "resolved_candidate"
    AMBIGUOUS_CANDIDATES = "ambiguous_candidates"
    INSUFFICIENT_COVERAGE = "insufficient_coverage"
    CONFLICTING_OBSERVATIONS = "conflicting_observations"
    NO_COMPATIBLE_DEFINITION = "no_compatible_definition"
    UNMODELED_OBSERVED_VARIATION = "unmodeled_observed_variation"
    UNSUPPORTED_GENE_OR_METHOD = "unsupported_gene_or_method"
    COMBINATORIAL_LIMIT_EXCEEDED = "combinatorial_limit_exceeded"


class PhenotypeStatus(StrEnum):
    """Source-rule mapping status; never a clinical recommendation."""

    MAPPED_UNIQUE_CANDIDATE = "mapped_unique_candidate"
    CONSISTENT_ACROSS_AMBIGUOUS_CANDIDATES = "consistent_across_ambiguous_candidates"
    DISCORDANT_ACROSS_CANDIDATES = "discordant_across_candidates"
    MAPPING_MISSING = "mapping_missing"
    SOURCE_VERSION_INCOMPATIBLE = "source_version_incompatible"
    NOT_EVALUABLE = "not_evaluable"


class CanonicalLocus(BaseModel):
    """Exact GRCh38 biallelic small-variant identity."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    locus_id: str
    gene_id: str
    assembly: str
    chromosome: str
    position: int = Field(gt=0)
    reference: str
    alternate: str


class LocusEvidence(BaseModel):
    """Observed evidence and lineage for one exact canonical locus."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    locus_id: str
    state: LocusEvidenceState
    alleles: tuple[str, ...] = ()
    genotype_ids: tuple[str, ...] = ()
    observation_references: tuple[str, ...] = ()


class CandidateDiplotype(BaseModel):
    """Canonical unordered software candidate, not a clinical genotype."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    candidate_id: str
    gene_id: str
    allele_a: str
    allele_b: str
    fully_evaluated: bool


class GeneInference(BaseModel):
    """Complete deterministic candidate set and conservative gene outcome."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    gene_id: str
    outcome: GeneOutcome
    candidates: tuple[CandidateDiplotype, ...] = ()
    locus_evidence: tuple[LocusEvidence, ...] = ()
    exclusion_reasons: tuple[str, ...] = ()


class MatchingLimits(BaseModel):
    """Validated bounds for the approximately O(A²L) generic matcher."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    max_alleles: int = Field(default=128, ge=1, le=4096)
    max_loci: int = Field(default=512, ge=1, le=100_000)
    max_candidate_pairs: int = Field(default=8256, ge=1, le=1_000_000)


class BundleValidation(BaseModel):
    """Validated immutable local bundle identity."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    directory: Path
    bundle_id: str
    bundle_version: str
    bundle_hash: str
    genes: tuple[str, ...]


class PharmacogenomicsResult(BaseModel):
    """Aggregate identity of an immutable M8 evidence run."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    run_id: str
    output_directory: Path
    gene_results: tuple[GeneInference, ...]
