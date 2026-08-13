from pathlib import Path

from pydantic import BaseModel, ConfigDict


class Settings(BaseModel):
    """Non-sensitive local settings. Data paths should normally be outside the repository."""

    model_config = ConfigDict(frozen=True)
    artifact_directory: Path | None = None
