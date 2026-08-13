"""Typed boundaries for transparent, non-classifying M4 review routing."""

from enum import StrEnum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class AnalysisContext(StrEnum):
    GERMLINE_CONSTITUTIONAL = "germline_constitutional"


class GenotypeEvidenceState(StrEnum):
    SINGLE_OBSERVED_ALT_CALL = "single_observed_alt_call"
    CONCORDANT_OBSERVED_ALT_CALLS = "concordant_observed_alt_calls"
    DISCORDANT_CALLED_ROWS_WITH_ALT = "discordant_called_rows_with_alt"
    OBSERVED_REFERENCE_ONLY = "observed_reference_only"
    NO_CANONICAL_CALLED_GENOTYPE = "no_canonical_called_genotype"


class CandidateEligibility(StrEnum):
    ELIGIBLE = "eligible"
    DATA_CONFLICT = "data_conflict"
    NOT_ELIGIBLE = "not_eligible"


class SourceTermRoute(StrEnum):
    HIGH_ATTENTION = "high_attention_source_term"
    RISK_CONTEXT = "risk_context"
    UNCERTAIN = "uncertain"
    BENIGN_LIKE = "benign_like"
    OTHER_CONTEXT = "other_context"
    UNMAPPED = "unmapped_source_term"


class ReviewPriorityBand(StrEnum):
    REVIEW_FIRST = "review_first"
    REVIEW_NEXT = "review_next"
    CONTEXT_ONLY = "context_only"
    NOT_PRIORITIZED = "not_prioritized"
    NOT_ELIGIBLE = "not_eligible"


class SourceReviewLevel(StrEnum):
    PRACTICE_GUIDELINE = "practice_guideline"
    EXPERT_PANEL = "expert_panel"
    MULTIPLE_SUBMITTERS = "multiple_submitters"
    SINGLE_SUBMITTER = "single_submitter"
    NO_ASSERTION_CRITERIA = "no_assertion_criteria"
    UNKNOWN = "unknown"


class PriorityRule(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    rule_id: str
    version: str
    band: ReviewPriorityBand


class PrioritizationPolicy(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    schema_version: int = Field(ge=1, le=1)
    policy_id: str
    version: str
    supported_analysis_context: AnalysisContext
    source_namespace: str
    dataset: str
    classification_type_handling: dict[str, str]
    term_routes: dict[SourceTermRoute, tuple[str, ...]]
    source_review_status_mappings: dict[str, SourceReviewLevel]
    record_status_behavior: dict[str, str]
    priority_rules: tuple[PriorityRule, ...]
    freshness_warning_days: int | None = Field(default=None, ge=0)
    tie_breaking: tuple[str, ...]

    @field_validator("version")
    @classmethod
    def semantic_version(cls, value: str) -> str:
        parts = value.split(".")
        if len(parts) != 3 or any(not part.isdigit() for part in parts):
            raise ValueError("policy version must be a three-part semantic version")
        return value

    @model_validator(mode="after")
    def routes_do_not_overlap(self) -> "PrioritizationPolicy":
        if set(self.term_routes) != set(SourceTermRoute):
            raise ValueError("policy must define every source-term route")
        expected_types = {"germline", "somatic_clinical_impact", "oncogenicity", "unknown"}
        if set(self.classification_type_handling) != expected_types:
            raise ValueError("policy must define every supported classification type")
        if self.classification_type_handling != {
            "germline": "route",
            "somatic_clinical_impact": "context_only",
            "oncogenicity": "context_only",
            "unknown": "context_only",
        }:
            raise ValueError("unsupported classification-type handling")
        if not self.priority_rules or len({r.rule_id for r in self.priority_rules}) != len(
            self.priority_rules
        ):
            raise ValueError("policy priority rule IDs must be present and unique")
        if any(
            value not in {"active", "inactive"} for value in self.record_status_behavior.values()
        ):
            raise ValueError("record-status behavior must be active or inactive")
        if set(self.record_status_behavior) != {"current", "replaced", "deleted", "unknown"}:
            raise ValueError("policy must define current, replaced, deleted, and unknown statuses")
        if (
            self.record_status_behavior["replaced"] != "inactive"
            or self.record_status_behavior["deleted"] != "inactive"
        ):
            raise ValueError("deleted and replaced assertions must be inactive")
        seen: dict[str, SourceTermRoute] = {}
        for route, terms in self.term_routes.items():
            for term in terms:
                normalized = " ".join(term.casefold().split())
                if normalized in seen:
                    raise ValueError(f"source term appears in multiple routes: {term}")
                seen[normalized] = route
        if self.tie_breaking[-1:] != ("profile_id",):
            raise ValueError("profile_id must be the final deterministic tie-breaker")
        return self


class PolicyIdentity(BaseModel):
    model_config = ConfigDict(frozen=True)
    policy_id: str
    version: str
    file_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    file_size_bytes: int = Field(ge=0)
    configuration_hash: str = Field(pattern=r"^[0-9a-f]{64}$")
    parsed_configuration: dict[str, Any]


class GenotypeRowEvidence(BaseModel):
    model_config = ConfigDict(frozen=True)
    genotype_id: str
    observation_reference: str
    ploidy: int
    alt_allele_count: int


class VariantEvidenceProfile(BaseModel):
    model_config = ConfigDict(frozen=True)
    profile_id: str
    m2_run_id: str
    variant_id: str
    assembly: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    genotype_state: GenotypeEvidenceState
    genotype_rows: tuple[GenotypeRowEvidence, ...]
    representation_ids: tuple[str, ...]
    assertion_ids: tuple[str, ...]
    scv_assertion_ids: tuple[str, ...]
    vcv_assertion_ids: tuple[str, ...]
    source_accessions: tuple[str, ...]
    assertion_levels: tuple[str, ...]
    classification_types: tuple[str, ...]
    source_terms: tuple[str, ...]
    source_term_routes: tuple[SourceTermRoute, ...]
    source_review_statuses: tuple[str, ...]
    source_review_levels: tuple[SourceReviewLevel, ...]
    source_record_statuses: tuple[str, ...]
    submitters: tuple[str, ...]
    condition_ids: tuple[str, ...]
    condition_names: tuple[str, ...]
    date_last_evaluated: tuple[str, ...]
    evidence_age_days: tuple[int | None, ...]
    source_snapshot_id: str
    source_release_date: str
    source_reported_conflict: bool
    submission_term_diversity: tuple[str, ...]
    missing_indicators: tuple[str, ...]
    unresolved_assessments: tuple[str, ...]


class ClinicalReviewCandidate(BaseModel):
    model_config = ConfigDict(frozen=True)
    candidate_id: str
    profile_id: str
    variant_id: str
    priority_band: ReviewPriorityBand
    eligibility: CandidateEligibility
    ordering_components: tuple[str, ...]
    source_review_level_order: int = Field(ge=0)
    source_reported_conflict_order: int = Field(ge=0, le=1)
    evaluation_date_missing_order: int = Field(ge=0, le=1)
    evidence_age_days_order: int = Field(ge=-1)


class CandidateAssertionLink(BaseModel):
    model_config = ConfigDict(frozen=True)
    link_id: str
    candidate_id: str
    assertion_instance_id: str
    assertion_level: str
    classification_type: str
    source_term: str
    source_term_route: SourceTermRoute
    active_for_routing: bool


class PriorityRationale(BaseModel):
    model_config = ConfigDict(frozen=True)
    rationale_id: str
    candidate_id: str
    profile_id: str
    policy_rule_id: str
    policy_rule_version: str
    rationale_type: str
    reason_code: str
    assertion_ids: tuple[str, ...] = ()
    genotype_ids: tuple[str, ...] = ()
    observation_references: tuple[str, ...] = ()
    source_term: str | None = None
    classification_type: str | None = None
    explanation: str


class PrioritizationExclusion(BaseModel):
    model_config = ConfigDict(frozen=True)
    exclusion_id: str
    profile_id: str
    candidate_id: str
    reason_code: str


class ClinicalPrioritizationConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    policy_path: Path
    analysis_context: AnalysisContext


class PrioritizationResult(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    run_id: str
    output_directory: Path
    policy_identity: PolicyIdentity
    profiles: tuple[VariantEvidenceProfile, ...]
    candidates: tuple[ClinicalReviewCandidate, ...]
    rationales: tuple[PriorityRationale, ...]
    candidate_assertion_links: tuple[CandidateAssertionLink, ...]
    exclusions: tuple[PrioritizationExclusion, ...]
    manifest: dict[str, Any]
