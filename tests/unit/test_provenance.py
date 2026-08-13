from datetime import UTC, datetime, timedelta

import pytest
from pydantic import ValidationError

from genome_evidence.domain.provenance import RunProvenance


def test_provenance_validates() -> None:
    started = datetime(2026, 1, 1, tzinfo=UTC)
    provenance = RunProvenance(
        run_id="run-1",
        started_at=started,
        completed_at=started + timedelta(seconds=1),
        git_commit="abc123",
        software_name="genome-evidence",
        software_version="0.1.0",
        configuration_hash="sha256:config",
        input_hashes={"fixture": "sha256:input"},
        transformation="synthetic test",
    )
    assert provenance.input_hashes["fixture"] == "sha256:input"


def test_provenance_rejects_reversed_timestamps() -> None:
    started = datetime(2026, 1, 2, tzinfo=UTC)
    with pytest.raises(ValidationError):
        RunProvenance(
            run_id="run-1",
            started_at=started,
            completed_at=started - timedelta(seconds=1),
            software_name="genome-evidence",
            software_version="0.1.0",
            configuration_hash="sha256:config",
            transformation="synthetic test",
        )
