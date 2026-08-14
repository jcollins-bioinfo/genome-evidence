import pytest

from genome_evidence.notebook_support import resolve_settings


def test_profile_resolution_is_explicit_and_safe() -> None:
    assert resolve_settings({}, in_colab=False).profile == "synthetic_ci"
    assert resolve_settings({}, in_colab=True).profile == "personal_drive"
    assert (
        resolve_settings(
            {"GENOME_EVIDENCE_PROFILE": "personal_drive", "GENOME_EVIDENCE_WORKSPACE": "/private"},
            in_colab=False,
        ).workspace_root.as_posix()
        == "/private"
    )
    with pytest.raises(ValueError):
        resolve_settings({"GENOME_EVIDENCE_PROFILE": "fallback"})
