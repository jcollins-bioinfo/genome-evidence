# Genome Evidence

Genome Evidence is an early-stage, reproducible personal-genomics **evidence infrastructure** project. It will preserve measurements, normalize representations, and integrate explicitly typed inferences and versioned evidence with end-to-end provenance. It is not a diagnostic service, clinical decision system, or a source of treatment recommendations.

## Epistemic boundary

An assay observation is not an inference. Measured genotypes, canonical alleles, imputed genotypes, external assertions, phenotype observations, and model outputs remain distinct. Missing means unknown or unassayed—never homozygous reference. Assay confidence and interpretation confidence are independent. Original observations are immutable; every transformation creates a traceable derived record. See the [epistemic contract](docs/epistemic-contract.md).

23andMe is planned as the first ingestion adapter, not the conceptual model. A later VCF/BCF/WGS adapter will enter through the same observation boundary, as will other assays, family samples, phenotypes, and laboratory measurements.

## Status

Milestone **M0 (bootstrap)** only: typed domain boundaries, storage schema scaffolding, documentation, and a local-only CLI. No ingestion, normalization, interpretation, imputation, scoring, reference downloads, or biological analysis is implemented.

## Install and use

Requires Python 3.12+ and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --dev
uv run genome-evidence version
uv run genome-evidence doctor
```

`doctor` checks only the local runtime and configuration; it does not inspect genotype data.

## Quality gates

```bash
uv run ruff check .
uv run ruff format --check .
uv run mypy src
uv run pytest
```

## Privacy warning

Never commit personal genomic, phenotype, clinical, or derived report data. Only clearly synthetic fixtures belong in tests. Prefer private data locations outside this repository. See [privacy guidance](docs/privacy.md). This project makes no regulatory-compliance claim.

## Roadmap

The staged plan runs from source ingestion and canonical normalization through versioned evidence, family and phenotype integration, and eventually WGS ingestion. Later milestones must not weaken the epistemic contract. See the [roadmap](docs/roadmap.md) and [architecture](docs/architecture.md).
