from importlib.metadata import version

from typer.testing import CliRunner

from genome_evidence.cli import app

runner = CliRunner()


def test_version() -> None:
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.stdout.strip() == version("genome-evidence")


def test_doctor_is_local_smoke_check() -> None:
    result = runner.invoke(app, ["doctor"])
    assert result.exit_code == 0
    assert f"Genome Evidence {version('genome-evidence')}: ok" in result.stdout
