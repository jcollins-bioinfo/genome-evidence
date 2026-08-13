"""Reusable provenance models."""

from datetime import datetime

from pydantic import BaseModel, ConfigDict, Field, model_validator


class ReferenceSourceVersion(BaseModel):
    model_config = ConfigDict(frozen=True)
    source: str = Field(min_length=1)
    version: str = Field(min_length=1)
    retrieved_at: datetime


class RunProvenance(BaseModel):
    model_config = ConfigDict(frozen=True)
    run_id: str = Field(min_length=1)
    started_at: datetime
    completed_at: datetime | None = None
    git_commit: str | None = None
    software_name: str = Field(min_length=1)
    software_version: str = Field(min_length=1)
    configuration_hash: str = Field(min_length=1)
    input_record_ids: tuple[str, ...] = ()
    input_hashes: dict[str, str] = Field(default_factory=dict)
    reference_sources: tuple[ReferenceSourceVersion, ...] = ()
    transformation: str = Field(min_length=1)

    @model_validator(mode="after")
    def completion_follows_start(self) -> "RunProvenance":
        if self.completed_at is not None and self.completed_at < self.started_at:
            raise ValueError("completed_at cannot precede started_at")
        return self
