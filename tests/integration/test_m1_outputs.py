import json
import subprocess
from pathlib import Path

import polars as pl
from typer.testing import CliRunner

from genome_evidence.cli import app
from genome_evidence.ingest import ingest_23andme

FIXTURE = Path(__file__).parents[1] / "fixtures" / "23andme" / "normal.txt"


def test_artifacts_and_manifest_are_provenance_linked(tmp_path: Path) -> None:
    output = tmp_path / "run"
    result = ingest_23andme(FIXTURE, output)
    expected = {
        "manifest.json",
        "source_metadata.json",
        "observations.parquet",
        "qc_summary.json",
        "qc_findings.parquet",
        "qc_report.md",
    }
    assert {p.name for p in output.iterdir()} == expected
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["run_id"] == result.run_id
    assert manifest["input_sha256"] == result.source_metadata.source_sha256
    assert manifest["configuration_hash"]
    assert set(manifest["artifacts"]) == expected - {"manifest.json"}
    assert pl.read_parquet(output / "observations.parquet")["raw_genotype"].to_list() == [
        "AG",
        "CC",
        "--",
    ]
    report = (output / "qc_report.md").read_text()
    assert "does not provide medical or biological interpretation" in report
    assert "rsSynthetic1" not in report


def test_substantive_artifact_hashes_repeat(tmp_path: Path) -> None:
    first = ingest_23andme(FIXTURE, tmp_path / "first")
    second = ingest_23andme(FIXTURE, tmp_path / "second")
    # Aggregate scientific content is stable; observation artifacts intentionally contain run IDs.
    for artifact in ("qc_summary.json", "qc_report.md"):
        assert first.manifest["artifacts"][artifact] == second.manifest["artifacts"][artifact]


def test_cli_ingestion_smoke(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        app, ["ingest", "23andme", "--input", str(FIXTURE), "--output", str(tmp_path / "run")]
    )
    assert result.exit_code == 0, result.output
    assert "3 source records" in result.output


def test_private_repository_output_pattern_is_ignored() -> None:
    result = subprocess.run(
        ["git", "check-ignore", "-q", "runs/private/synthetic-run/manifest.json"], check=False
    )
    assert result.returncode == 0
