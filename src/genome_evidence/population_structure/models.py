"""Typed public M5 population-structure contracts."""

from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, model_validator


class FrozenModel(BaseModel):
    model_config = ConfigDict(frozen=True)


class BundleArtifact(FrozenModel):
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    byte_size: int = Field(ge=0)


class PopulationReferenceIdentity(FrozenModel):
    model_id: str
    model_version: str
    assembly: str
    manifest_sha256: str
    artifacts: dict[str, BundleArtifact]


class ReferenceBundleManifest(FrozenModel):
    schema_version: int
    model_id: str
    model_version: str
    method: str
    algorithm_version: str
    assembly: str
    supported_chromosomes: tuple[str, ...]
    marker_count: int = Field(gt=0)
    reference_sample_count: int = Field(gt=0)
    reference_group_count: int = Field(gt=0)
    component_count: int = Field(gt=0)
    default_component_count: int = Field(gt=0)
    maximum_component_count: int = Field(gt=0)
    source_datasets: tuple[str, ...]
    source_releases: tuple[str, ...]
    citations: tuple[str, ...]
    stable_urls: tuple[str, ...]
    creation_date: str
    licensing: str
    redistribution_terms: str
    reference_label_definitions: str
    reference_label_provenance: str
    training_methodology: str
    preprocessing_methodology: str
    variant_qc: str
    missingness_rules: str
    allele_frequency_handling: str
    ld_pruning: str
    relatedness_handling: str
    sample_qc_and_outliers: str
    sample_imbalance_handling: str
    standardization_convention: str
    pca_implementation: str
    matrix_convention: str
    axis_sign_convention: str
    projection_method: str
    least_squares_rcond: float = Field(gt=0)
    numerical_tolerance: float = Field(gt=0, le=1e-6)
    minimum_observed_marker_count: int = Field(gt=0)
    minimum_observed_marker_fraction: float = Field(gt=0, le=1)
    minimum_chromosomes: int = Field(gt=0, le=22)
    minimum_loading_energy: float = Field(gt=0, le=1)
    maximum_condition_number: float = Field(ge=1)
    support_envelope_method: str
    support_quantile: float = Field(gt=0, lt=1)
    support_quantile_method: str
    support_minimum_group_size: int = Field(gt=1)
    nearest_neighbor_distance_convention: str
    sensitivity_minimum_valid_replicates: int = Field(gt=0)
    projection_validation_scope: str
    known_limitations: str
    intended_use: str
    explicit_non_goals: str
    artifacts: dict[str, BundleArtifact]

    @model_validator(mode="after")
    def conventions(self) -> "ReferenceBundleManifest":
        if self.schema_version != 1 or self.assembly != "GRCh38":
            raise ValueError("unsupported reference schema or assembly")
        if self.matrix_convention != "orthonormal_marker_by_component_float64":
            raise ValueError("unsupported loading matrix convention")
        if self.projection_method != "partial_marker_least_squares":
            raise ValueError("unsupported projection method")
        if self.default_component_count > self.maximum_component_count > self.component_count:
            raise ValueError("invalid component limits")
        required = {
            "variant_loadings.parquet",
            "reference_scores.parquet",
            "reference_groups.parquet",
            "component_metadata.parquet",
        }
        if set(self.artifacts) != required:
            raise ValueError("reference artifact declaration is incomplete")
        return self


class ReferenceVariantLoading(FrozenModel):
    model_variant_id: str
    assembly: str
    chromosome: str
    position: int = Field(gt=0)
    reference: str
    alternate: str
    effect_allele: str
    training_alt_allele_frequency: float = Field(ge=0, le=1)
    training_mean_dosage: float = Field(ge=0, le=2)
    training_scale: float = Field(gt=0)
    loadings: tuple[float, ...]


class ReferenceGroup(FrozenModel):
    reference_group_id: str
    source_label: str
    source_definition: str
    label_provenance: str
    citation: str
    declared_sample_count: int = Field(gt=0)
    parent_group: str | None = None


class ComponentMetadata(FrozenModel):
    component_index: int = Field(gt=0)
    eigenvalue: float = Field(gt=0)
    explained_variance_fraction: float = Field(gt=0, le=1)
    reference_score_mean: float
    reference_score_standard_deviation: float = Field(gt=0)
    sign_calibration: str


class ReferenceValidationResult(FrozenModel):
    directory: Path
    identity: PopulationReferenceIdentity
    manifest: ReferenceBundleManifest
    variants: tuple[ReferenceVariantLoading, ...]
    groups: tuple[ReferenceGroup, ...]
    components: tuple[ComponentMetadata, ...]
    reference_scores: tuple[tuple[str, str, tuple[float, ...]], ...]


class PopulationStructureConfig(FrozenModel):
    component_count: int | None = Field(default=None, gt=0)
    nearest_neighbor_count: int = Field(default=10, gt=0)
    minimum_observed_marker_count: int | None = Field(default=None, gt=0)
    minimum_observed_marker_fraction: float | None = Field(default=None, gt=0, le=1)
    minimum_chromosomes: int | None = Field(default=None, gt=0, le=22)
    minimum_loading_energy: float | None = Field(default=None, gt=0, le=1)
    maximum_condition_number: float | None = Field(default=None, ge=1)
    sensitivity_minimum_valid_replicates: int | None = Field(default=None, gt=0)


class MarkerAlignmentOutcome(StrEnum):
    USED = "used"
    MISSING = "missing"
    EXCLUDED = "excluded"
    DISCORDANT = "discordant"


class MarkerExclusionReason(StrEnum):
    MODEL_MARKER_ABSENT = "model_marker_absent"
    DISCORDANT_CALLS = "discordant_calls"
    NONAUTOSOMAL = "nonautosomal"
    UNSUPPORTED_VARIANT = "unsupported_variant"
    HAPLOID = "haploid"
    INVALID_ALLELES = "invalid_alleles"


class MarkerAlignment(FrozenModel):
    alignment_id: str
    model_variant_id: str
    variant_id: str | None
    outcome: MarkerAlignmentOutcome
    exclusion_reason: MarkerExclusionReason | None
    alt_dosage: int | None = Field(default=None, ge=0, le=2)
    genotype_ids: tuple[str, ...] = ()
    observation_references: tuple[str, ...] = ()


class ProjectionStatus(StrEnum):
    PROJECTED = "projected"
    NOT_PROJECTED = "not_projected"


class ProjectionDiagnostic(FrozenModel):
    selected_component_count: int
    used_marker_count: int
    observed_marker_fraction: float
    chromosomes: tuple[str, ...]
    loading_energy_coverage: tuple[float, ...]
    matrix_rank: int
    singular_values: tuple[float, ...]
    condition_number: float | None
    residual_norm: float | None
    failed_gates: tuple[str, ...]


class ReferenceSupportStatus(StrEnum):
    WITHIN = "within_at_least_one_reference_envelope"
    OUTSIDE = "outside_all_evaluable_reference_envelopes"
    NOT_EVALUABLE = "not_evaluable"


class ProjectionCoordinate(FrozenModel):
    coordinate_id: str
    component_index: int
    coordinate: float


class ReferenceGroupDistance(FrozenModel):
    reference_group_id: str
    source_label: str
    source_definition: str
    label_provenance: str
    sample_count: int
    distance: float
    rank: int


class ReferenceNeighbor(FrozenModel):
    reference_sample_id: str
    reference_group_id: str
    distance: float
    rank: int


class ReferenceSupportEvaluation(FrozenModel):
    reference_group_id: str
    evaluable: bool
    inside_envelope: bool | None
    distance: float | None
    envelope_radius: float | None


class ProjectionSensitivityReplicate(FrozenModel):
    omitted_chromosome: str
    status: ProjectionStatus
    retained_marker_count: int
    retained_marker_fraction: float
    diagnostic: ProjectionDiagnostic
    coordinates: tuple[float, ...] = ()
    absolute_deviations: tuple[float, ...] = ()
    rank_one_group_id: str | None = None
    support_status: ReferenceSupportStatus | None = None


class PopulationStructureResult(FrozenModel):
    run_id: str
    output_directory: Path
    reference_identity: PopulationReferenceIdentity
    projection_status: ProjectionStatus
    support_status: ReferenceSupportStatus | None
    alignments: tuple[MarkerAlignment, ...]
    diagnostic: ProjectionDiagnostic
    coordinates: tuple[ProjectionCoordinate, ...]
    group_distances: tuple[ReferenceGroupDistance, ...]
    neighbors: tuple[ReferenceNeighbor, ...]
    support: tuple[ReferenceSupportEvaluation, ...]
    sensitivity: tuple[ProjectionSensitivityReplicate, ...]
