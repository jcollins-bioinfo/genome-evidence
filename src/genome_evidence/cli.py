import platform
from pathlib import Path
from typing import Annotated

import typer

from genome_evidence import __version__
from genome_evidence.evidence import (
    AnnotationConfig,
    ClinVarIngestionConfig,
    ingest_clinvar_vcv,
    link_external_evidence,
)
from genome_evidence.ingest import Ingest23andMeConfig, ParseMode, ingest_23andme
from genome_evidence.ingest.errors import GenotypeParseError
from genome_evidence.normalization import NormalizationConfig, normalize_m1_run
from genome_evidence.population_structure import (
    PopulationStructureConfig,
    infer_population_structure,
    validate_population_reference,
)
from genome_evidence.prioritization import prioritize_clinical_variants
from genome_evidence.prioritization.models import AnalysisContext, ClinicalPrioritizationConfig

app = typer.Typer(help="Genome Evidence utilities.", no_args_is_help=True)
ingest_app = typer.Typer(help="Ingest source-faithful observations.", no_args_is_help=True)
app.add_typer(ingest_app, name="ingest")
evidence_app = typer.Typer(
    help="Ingest and link versioned external assertions.", no_args_is_help=True
)
app.add_typer(evidence_app, name="evidence")
prioritize_app = typer.Typer(help="Create transparent manual-review queues.", no_args_is_help=True)
app.add_typer(prioritize_app, name="prioritize")
ancestry_app = typer.Typer(
    help="Reference-panel population-structure projection.", no_args_is_help=True
)
app.add_typer(ancestry_app, name="ancestry")


@ancestry_app.command("validate-reference")
def ancestry_validate_reference(
    reference_bundle: Annotated[
        Path, typer.Option("--reference-bundle", exists=True, file_okay=False)
    ],
) -> None:
    """Validate one local, checksummed population reference bundle."""
    try:
        result = validate_population_reference(reference_bundle)
    except (ValueError, OSError) as error:
        typer.echo(f"Reference validation failed: {error}", err=True)
        raise typer.Exit(code=2) from error
    typer.echo(f"Reference valid: {result.identity.model_id} {result.identity.model_version}")
    typer.echo(
        f"Markers: {len(result.variants)}; reference samples: {len(result.reference_scores)}"
    )


@ancestry_app.command("project")
def ancestry_project(
    normalization_run: Annotated[
        Path, typer.Option("--normalization-run", exists=True, file_okay=False)
    ],
    reference_bundle: Annotated[
        Path, typer.Option("--reference-bundle", exists=True, file_okay=False)
    ],
    output: Annotated[Path, typer.Option("--output", file_okay=False)],
    components: Annotated[int | None, typer.Option("--components", min=1)] = None,
    neighbors: Annotated[int, typer.Option("--neighbors", min=1)] = 10,
) -> None:
    """Project one validated M2 analysis unit into a fixed local PCA space."""
    try:
        result = infer_population_structure(
            normalization_run,
            reference_bundle,
            output,
            PopulationStructureConfig(component_count=components, nearest_neighbor_count=neighbors),
        )
    except (ValueError, FileExistsError, OSError) as error:
        typer.echo(f"Population-structure projection failed: {error}", err=True)
        raise typer.Exit(code=2) from error
    typer.echo(f"Projection status: {result.projection_status}")
    typer.echo(f"Support status: {result.support_status}")
    typer.echo(f"Used markers: {result.diagnostic.used_marker_count}")
    typer.echo(f"Output: {result.output_directory}")
    typer.echo(f"Run ID: {result.run_id}")
    typer.echo(
        f"Reference: {result.reference_identity.model_id} {result.reference_identity.model_version}"
    )


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


@evidence_app.command("ingest-clinvar")
def evidence_ingest_clinvar(
    input_path: Annotated[Path, typer.Option("--input", exists=True, dir_okay=False)],
    output_path: Annotated[Path, typer.Option("--output", file_okay=False)],
    release_identity: Annotated[str | None, typer.Option("--release-identity")] = None,
) -> None:
    """Ingest one fixed, local ClinVar VCV XML or XML.GZ snapshot."""
    try:
        result = ingest_clinvar_vcv(
            input_path,
            output_path,
            ClinVarIngestionConfig(release_identity_override=release_identity),
        )
    except (ValueError, FileExistsError, OSError) as error:
        typer.echo(f"Evidence ingestion failed: {error}", err=True)
        raise typer.Exit(code=2) from error
    typer.echo(f"Evidence ingestion complete: {len(result.assertions)} assertions")
    typer.echo(f"Run ID: {result.run_id}")


@evidence_app.command("link")
def evidence_link(
    normalization_run: Annotated[
        Path, typer.Option("--normalization-run", exists=True, file_okay=False)
    ],
    evidence_run: Annotated[Path, typer.Option("--evidence-run", exists=True, file_okay=False)],
    output_path: Annotated[Path, typer.Option("--output", file_okay=False)],
) -> None:
    """Link external representations to M2 by exact canonical allele identity."""
    try:
        result = link_external_evidence(
            normalization_run, evidence_run, output_path, AnnotationConfig()
        )
    except (ValueError, FileExistsError, OSError) as error:
        typer.echo(f"Evidence linking failed: {error}", err=True)
        raise typer.Exit(code=2) from error
    matched = sum(link.outcome.value == "matched" for link in result.links)
    typer.echo(f"Evidence linking complete: {len(result.links)} representations; {matched} matched")
    typer.echo(f"Run ID: {result.run_id}")


@prioritize_app.command("clinical")
def prioritize_clinical(
    normalization_run: Annotated[
        Path, typer.Option("--normalization-run", exists=True, file_okay=False)
    ],
    evidence_run: Annotated[Path, typer.Option("--evidence-run", exists=True, file_okay=False)],
    annotation_run: Annotated[Path, typer.Option("--annotation-run", exists=True, file_okay=False)],
    policy: Annotated[Path, typer.Option("--policy", exists=True, dir_okay=False)],
    analysis_context: Annotated[AnalysisContext, typer.Option("--analysis-context")],
    output: Annotated[Path, typer.Option("--output", file_okay=False)],
) -> None:
    """Route source-linked records for manual review using only local inputs."""
    try:
        result = prioritize_clinical_variants(
            normalization_run,
            evidence_run,
            annotation_run,
            output,
            ClinicalPrioritizationConfig(policy_path=policy, analysis_context=analysis_context),
        )
    except (ValueError, FileExistsError, OSError) as error:
        typer.echo(f"Prioritization failed: {error}", err=True)
        raise typer.Exit(code=2) from error
    counts: dict[str, int] = {}
    for candidate in result.candidates:
        counts[candidate.priority_band.value] = counts.get(candidate.priority_band.value, 0) + 1
    typer.echo(
        f"Prioritization complete: {len(result.profiles)} profiles; "
        f"{len(result.candidates)} candidates"
    )
    typer.echo("Bands: " + ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    typer.echo(f"Output: {result.output_directory}")
    typer.echo(f"Run ID: {result.run_id}")


if __name__ == "__main__":
    app()
