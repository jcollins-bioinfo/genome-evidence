from typer.testing import CliRunner

from genome_evidence.cli import app

runner = CliRunner()


def test_version() -> None:
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.stdout.strip() == "0.1.0"


def test_doctor_is_local_smoke_check() -> None:
    result = runner.invoke(app, ["doctor"])
    assert result.exit_code == 0
    assert "Genome Evidence 0.1.0: ok" in result.stdout
