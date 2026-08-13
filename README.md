# Genome Evidence

Genome Evidence is an early-stage, reproducible personal-genomics **evidence infrastructure** project. It will preserve measurements, normalize representations, and integrate explicitly typed inferences and versioned evidence with end-to-end provenance. It is not a diagnostic service, clinical decision system, or a source of treatment recommendations.

## Epistemic boundary

An assay observation is not an inference. Measured genotypes, canonical alleles, imputed genotypes, external assertions, phenotype observations, and model outputs remain distinct. Missing means unknown or unassayed—never homozygous reference. Assay confidence and interpretation confidence are independent. Original observations are immutable; every transformation creates a traceable derived record. See the [epistemic contract](docs/epistemic-contract.md).

23andMe is planned as the first ingestion adapter, not the conceptual model. A later VCF/BCF/WGS adapter will enter through the same observation boundary, as will other assays, family samples, phenotypes, and laboratory measurements.

## Status

Milestones **M1–M2** are implemented. M3 and M4 code are present; locked dependency and independent CI verification remain blocked until dependency resolution is available. M4 creates only a manual-review queue, never an automatic classification. M5 and later milestones remain unimplemented.

## Install and use

Requires Python 3.12+ and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --dev
uv run genome-evidence version
uv run genome-evidence doctor
uv run genome-evidence ingest 23andme --input /private/input.txt --output /private/m1-run
uv run genome-evidence normalize --input /private/m1-run --output /private/m2-run --marker-definitions /references/markers.json --target-reference /references/GRCh38.fa
uv run genome-evidence evidence ingest-clinvar --input /references/fixed-vcv.xml.gz --output /private/evidence-run
uv run genome-evidence evidence link --normalization-run /private/m2-run --evidence-run /private/evidence-run --output /private/annotation-run
uv run genome-evidence prioritize clinical --normalization-run /private/m2-run --evidence-run /private/evidence-run --annotation-run /private/annotation-run --policy references/clinvar-germline-review-policy-v1.json --analysis-context germline_constitutional --output /private/m4-run
```

`doctor` checks only the local runtime and configuration; it does not inspect genotype data. The ingestion CLI is a thin wrapper over `ingest_23andme`; see [the M1 ingestion guide](docs/ingestion/23andme.md) and [assay-QC definitions](docs/qc/assay-level-qc.md).

## Quality gates

```bash
uv run ruff check .
uv run ruff format --check .
uv run mypy src
uv run pytest
uv run pytest --nbmake notebooks
```

## Privacy warning

Never commit personal genomic, phenotype, clinical, or derived report data. Only clearly synthetic fixtures belong in tests. Prefer private data locations outside this repository. See [privacy guidance](docs/privacy.md). This project makes no regulatory-compliance claim.

## Roadmap

The staged plan runs from source ingestion and canonical normalization through versioned evidence, family and phenotype integration, and eventually WGS ingestion. Later milestones must not weaken the epistemic contract. See the [roadmap](docs/roadmap.md) and [architecture](docs/architecture.md).
