from datetime import datetime
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from genome_evidence.domain.subjects import Subject


class Sample(BaseModel):
    model_config = ConfigDict(frozen=True)
    sample_id: str = Field(min_length=1)
    subject: Subject
    source: str = Field(min_length=1)
    assay_type: str = Field(min_length=1)
    genome_build: str = Field(min_length=1)
    source_identifier: str = Field(min_length=1)
    ingested_at: datetime
    ingestion_metadata: dict[str, Any] = Field(default_factory=dict)
    source_checksum: str = Field(min_length=1)
