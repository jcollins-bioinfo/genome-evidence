"""Typed M9 contracts for declared pedigrees and site-level family evidence.

The models deliberately separate user-declared relationship assertions, directly
observed M2 genotype evidence, deterministic Mendelian compatibility, and
site-local transmission constraints.  None of these records verifies relatedness
or supports a clinical conclusion.
"""

from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, model_validator

PEDIGREE_SCHEMA = "genome-evidence-pedigree/v1"
RUN_SCHEMA = "genome-evidence-m9-family-run/v1"
ALGORITHM_VERSION = "m9-site-transmission-enumeration-v1"


class RelationshipSource(StrEnum):
    """Non-biological category describing how a relationship was declared."""

    USER_DECLARED = "user_declared"
    PROJECT_RECORD = "project_record"


class CompatibilityStatus(StrEnum):
    """Conditional Mendelian compatibility or a reason it was not evaluated."""

    CONSISTENT = "consistent"
    INCONSISTENT = "inconsistent"
    INDETERMINATE_MISSING = "indeterminate_missing"
    INDETERMINATE_CONFLICTING = "indeterminate_conflicting"
    UNSUPPORTED_LOCUS = "unsupported_locus"
    UNSUPPORTED_PLOIDY = "unsupported_ploidy"
    UNSUPPORTED_GENOTYPE_REPRESENTATION = "unsupported_genotype_representation"


class Informativeness(StrEnum):
    """Whether observations constrain a transmission beyond compatibility."""

    INFORMATIVE = "informative"
    UNINFORMATIVE = "uninformative"
    NOT_EVALUATED = "not_evaluated"


class TransmissionStatus(StrEnum):
    """Scope-limited status for transmissions at exactly one canonical site."""

    UNIQUE_TRANSMISSION = "unique_transmission"
    AMBIGUOUS_TRANSMISSION = "ambiguous_transmission"
    NO_COMPATIBLE_TRANSMISSION = "no_compatible_transmission"
    UNRESOLVED_INPUT = "unresolved_input"
    NOT_APPLICABLE = "not_applicable"


class PedigreeMember(BaseModel):
    """A pseudonymous member explicitly bound to one completed M2 run."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    member_id: str = Field(min_length=1, pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
    subject_id: str = Field(pattern=r"^subject-[0-9]{4,}$")
    m2_run: Path


class ParentRelationship(BaseModel):
    """A declared biological-parent assertion, not verified relatedness."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    assertion_id: str = Field(min_length=1)
    parent_member_id: str = Field(min_length=1)
    child_member_id: str = Field(min_length=1)
    source: RelationshipSource
    declared_role: str | None = None


class PedigreeDescriptor(BaseModel):
    """Strict local M9 input containing no names or clinical narrative."""

    model_config = ConfigDict(frozen=True, extra="forbid")
    schema_id: str = Field(default=PEDIGREE_SCHEMA, alias="schema")
    family_id: str = Field(min_length=1, pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
    members: tuple[PedigreeMember, ...]
    relationships: tuple[ParentRelationship, ...]

    @model_validator(mode="after")
    def known_schema(self) -> "PedigreeDescriptor":
        """Reject unknown schemas before any scientific evaluation."""
        if self.schema_id != PEDIGREE_SCHEMA:
            raise ValueError("unsupported pedigree schema")
        return self


class TransmissionAssignment(BaseModel):
    """One distinct site-local allele assignment in neutral parent-edge order."""

    model_config = ConfigDict(frozen=True)
    parent_1_allele: str
    parent_2_allele: str


class SegregationEvidence(BaseModel):
    """Deterministic evidence conditional on declared edges and observed calls."""

    model_config = ConfigDict(frozen=True)
    evidence_id: str
    variant_id: str
    child_member_id: str
    relationship_assertion_ids: tuple[str, ...]
    genotype_record_ids: tuple[str, ...]
    compatibility: CompatibilityStatus
    informativeness: Informativeness
    transmission_status: TransmissionStatus
    assignments: tuple[TransmissionAssignment, ...] = ()


class FamilyAnalysisResult(BaseModel):
    """Aggregate handle returned after immutable M9 artifact publication."""

    model_config = ConfigDict(frozen=True)
    run_id: str
    output_directory: Path
    evidence_count: int
    compatibility_counts: dict[str, int]
