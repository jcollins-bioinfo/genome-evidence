"""Package-version single-authority regression tests."""

import tomllib
from importlib.metadata import version
from pathlib import Path

import genome_evidence


def test_pyproject_distribution_and_runtime_versions_agree() -> None:
    project = tomllib.loads(Path("pyproject.toml").read_text())["project"]
    assert project["version"] == version("genome-evidence") == genome_evidence.__version__
