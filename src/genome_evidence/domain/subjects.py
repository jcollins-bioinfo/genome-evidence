from pydantic import BaseModel, ConfigDict, Field


class Subject(BaseModel):
    """A person represented only by a project-local pseudonym."""

    model_config = ConfigDict(frozen=True)
    subject_id: str = Field(min_length=1, description="Pseudonymous identifier")
