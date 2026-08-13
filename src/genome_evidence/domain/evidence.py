from datetime import date
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from genome_evidence.domain.enums import EvidenceType, InterpretationStatus


class EntityReference(BaseModel):
    model_config = ConfigDict(frozen=True)
    entity_type: str = Field(min_length=1)
    identifier: str = Field(min_length=1)


class EvidenceAssertion(BaseModel):
    """A versioned directed assertion; status is not a clinical grading system."""

    model_config = ConfigDict(frozen=True)
    assertion_id: str = Field(min_length=1)
    subject: EntityReference
    predicate: str = Field(min_length=1)
    object: EntityReference
    evidence_type: EvidenceType
    interpretation_status: InterpretationStatus
    source: str = Field(min_length=1)
    source_accession: str = Field(min_length=1)
    source_version: str = Field(min_length=1)
    retrieved_on: date
    evidence_metadata: dict[str, Any] = Field(default_factory=dict)
    run_id: str = Field(min_length=1)
