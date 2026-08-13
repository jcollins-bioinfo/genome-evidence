import platform
from pathlib import Path
from typing import Annotated

import typer

from genome_evidence import __version__
from genome_evidence.ingest import Ingest23andMeConfig, ParseMode, ingest_23andme
from genome_evidence.ingest.errors import GenotypeParseError
from genome_evidence.normalization import NormalizationConfig, normalize_m1_run

app = typer.Typer(help="Genome Evidence utilities.", no_args_is_help=True)
ingest_app = typer.Typer(help="Ingest source-faithful observations.", no_args_is_help=True)
app.add_typer(ingest_app, name="ingest")


@app.command()
def version() -> None:
    """Print the installed package version."""
    typer.echo(__version__)


@app.command()
def doctor() -> None:
    """Check only local runtime assumptions; never inspect genomic data."""
    current = (platform.python_version_tuple()[0], platform.python_version_tuple()[1])
    supported = tuple(map(int, current)) >= (3, 12)
    typer.echo(f"Python {platform.python_version()}: {'ok' if supported else 'unsupported'}")
    typer.echo(f"Genome Evidence {__version__}: ok")
    if not supported:
        raise typer.Exit(code=1)


@ingest_app.command("23andme")
def ingest_twenty_three_and_me(
    input_path: Annotated[
        Path, typer.Option("--input", exists=True, dir_okay=False, readable=True)
    ],
    output_path: Annotated[Path, typer.Option("--output", file_okay=False)],
    mode: Annotated[ParseMode, typer.Option()] = ParseMode.STRICT,
    genome_build: Annotated[str | None, typer.Option("--genome-build")] = None,
    sample_id: Annotated[str, typer.Option("--sample-id")] = "sample",
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Ingest a 23andMe raw file without biological normalization."""
    try:
        result = ingest_23andme(
            input_path,
            output_path,
            Ingest23andMeConfig(
                mode=mode,
                genome_build_override=genome_build,
                sample_id=sample_id,
                overwrite=overwrite,
            ),
        )
    except (GenotypeParseError, FileExistsError, NotADirectoryError, ValueError) as error:
        typer.echo(f"Ingestion failed: {error}", err=True)
        raise typer.Exit(code=2) from error
    typer.echo(f"Ingestion complete: {result.qc_summary.parsed_record_count} source records")
    typer.echo(f"Run ID: {result.run_id}")


@app.command("normalize")
def normalize(
    input_path: Annotated[Path, typer.Option("--input", exists=True, file_okay=False)],
    output_path: Annotated[Path, typer.Option("--output", file_okay=False)],
    marker_definitions: Annotated[Path, typer.Option("--marker-definitions", exists=True)],
    target_reference: Annotated[Path, typer.Option("--target-reference", exists=True)],
    target_build: Annotated[str, typer.Option("--target-build")] = "GRCh38",
    liftover: Annotated[Path | None, typer.Option("--liftover")] = None,
    source_build: Annotated[str | None, typer.Option("--source-build")] = None,
) -> None:
    """Normalize a completed M1 run using explicit local reference resources."""
    try:
        result = normalize_m1_run(
            input_path,
            output_path,
            NormalizationConfig(
                marker_definitions=marker_definitions,
                target_reference=target_reference,
                target_build=target_build,
                liftover=liftover,
                source_build_override=source_build,
            ),
        )
    except (ValueError, FileExistsError, OSError) as error:
        typer.echo(f"Normalization failed: {error}", err=True)
        raise typer.Exit(code=2) from error
    mapped = sum(mapping.outcome.value == "mapped" for mapping in result.mappings)
    typer.echo(f"Normalization complete: {len(result.mappings)} mappings; {mapped} mapped")
    typer.echo(f"Output: {result.output_directory}")
    typer.echo(f"Run ID: {result.run_id}")


if __name__ == "__main__":
    app()
