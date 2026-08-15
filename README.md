# Genome Evidence

Genome Evidence is early-stage, reproducible, provenance-first personal-genomics **evidence infrastructure**. It preserves observations, typed inferences, source-versioned assertions, uncertainty, and lineage. It is not a diagnostic service, clinical test, prescribing engine, or treatment recommendation system.

## Epistemic and privacy boundary

An observation is not an inference; a canonical allele is not a source observation; an external assertion is not project truth; a model candidate is not a clinical result. Missing means unknown or unassayed—never reference. Keep all personal genomic, phenotype, medication, clinical, and derived data outside this repository and offline where possible. Colab runs on a Google-managed VM. See the [epistemic contract](docs/epistemic-contract.md) and [privacy guidance](docs/privacy.md).

## Status through M8

M1 ingestion, M2 canonical normalization, M3 external evidence, M4 transparent review prioritization, M5 population projection foundations, M7 polygenic-score foundations, and M8 pharmacogenomic evidence foundations are implemented with synthetic tests. M6 remains gated: its infrastructure exists, but the required real-Beagle fabricated-panel smoke has not established completed real-engine support. No production population, imputation, PGS, or PGx reference bundle ships.

M8 validates strict local bundles and enumerates conservative candidates from exact M2 **observed** GRCh38 genotypes. Sparse consumer arrays are not complete pharmacogene assays. Missing positions never default to `*1`; ambiguity is retained. CYP2D6, CNV/SV/hybrids, HLA, mitochondrial, sex-chromosome, special-algorithm, and unsupported-ploidy calling fail closed. M8 produces no clinical validation, medication selection, dose, safety label, actionability claim, or recommendation.

## Install and quickstart

Requires Python 3.12+ and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --locked --dev
uv run genome-evidence version
uv run genome-evidence doctor
uv run genome-evidence --help
```

Selected command groups include `ingest`, `normalize`, `evidence`, `prioritize`, `ancestry`, `phasing`, `pgs`, `pgx`, and `workspace`. PGx analysis is offline and requires explicit genes:

```bash
uv run genome-evidence pgx validate-bundle --bundle /references/pgx-bundle
uv run genome-evidence pgx infer --normalization-run /private/m2-run --bundle /references/pgx-bundle --output /private/m8-run --gene FAKE1
```

These commands report aggregate status only. Use the public Python APIs `validate_pharmacogenomics_bundle(...)` and `infer_pharmacogenomics(...)` for the same checked workflow.

## Private workspace and artifacts

`genome-evidence workspace init` creates a provider-neutral, non-destructive tree. Durable PGx bundles belong under `references/pharmacogenomics/`, evictable source caches under `cache/{clinpgx,pharmvar,pharmcat}/`, and private M8 outputs under `runs/m8_pharmacogenomics/`. The canonical Colab location remains `/content/drive/MyDrive/genome-evidence-private`; local/offline roots are supported and preferred for privacy-sensitive work.

Notebook 00B provisions the production normalization resources required by the source imported in notebook 00. It downloads the pinned ~1.8 GiB `dbSnp155Common.bb` for each required assembly with resumable ranges, queries a verified local `/content` copy, and sends only common-missing or validation-indeterminate rsIDs through a bounded parallel authoritative query of the complete remote index. It never downloads the 65/68 GiB complete BigBeds and never transmits genotype calls. It verifies UCSC's GRCh38 FASTA checksum, builds the FAI and any conservative variant-specific GRCh37→GRCh38 map, and publishes provenance plus durable selectors last. Notebook 01 consumes that selection, computes ephemerally, and publishes re-hashed M1/M2 runs.

00B is restartable at operational boundaries. It writes privacy-safe progress both to the
notebook cell and to
`logs/notebooks/00b/<run-key>/events.jsonl`, while verified download and dbSNP-query
checkpoints live under `cache/downloads/normalization/v1/<run-key>/`. Rerunning the
provisioning cell verifies and reuses completed batches and downloaded bytes. A rejected common-query leaf is localized: validated sibling batches remain usable, while only that leaf and genuine common misses enter full-index fallback. Full-file
transfers and byte-processing stages report byte rates and ETA; sparse remote BigBed
queries report identifiers per second because the Kent utility does not expose truthful
per-request byte telemetry. Last-written bundle and selector completion markers make
pre-commit Drive/FUSE interruptions repairable without weakening post-commit
immutability. A retained `hg38.fa.gz` checkpoint adds about 938 MiB to the
approximately 3 GiB installed FASTA footprint, before smaller query checkpoints and logs.

Every stage creates immutable, checksummed artifacts with source/input identities, configuration and algorithm identity, package/Git identity, and explicit schemas. Package versions, artifact schemas, project milestones, and external source releases are independent. See [architecture](docs/architecture.md), the [data model](docs/data-model.md), [PGx model documentation](docs/pharmacogenomics/model-and-evidence.md), [changelog](CHANGELOG.md), and [version policy](docs/development/versioning-and-releases.md).

## Notebooks

The [canonical notebook index](notebooks/README.md) lists all nine notebooks exactly once. A badge only opens a notebook; it does not auto-run, mount Drive, grant access, create folders, or reproduce the lockfile environment.

- [00 initialize private workspace](notebooks/00_initialize_private_workspace.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/00_initialize_private_workspace.ipynb)
- [00B provision normalization resources](notebooks/00b_provision_normalization_resources.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/00b_provision_normalization_resources.ipynb)
- [01 ingest and normalize genome](notebooks/01_ingest_and_normalize_genome.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/01_ingest_and_normalize_genome.ipynb)
- [02 versioned external evidence](notebooks/02_versioned_external_evidence.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/02_versioned_external_evidence.ipynb)
- [03 evidence-oriented clinical review](notebooks/03_evidence_oriented_clinical_review.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/03_evidence_oriented_clinical_review.ipynb)
- [04 population structure projection](notebooks/04_population_structure_projection.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/04_population_structure_projection.ipynb)
- [05 reference-panel phasing and imputation](notebooks/05_reference_panel_phasing_and_imputation.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/05_reference_panel_phasing_and_imputation.ipynb)
- [06 versioned polygenic-score calculation](notebooks/06_polygenic_score_calculation.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/06_polygenic_score_calculation.ipynb)
- [07 pharmacogenomic candidate evidence](notebooks/07_pharmacogenomics_evidence.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/07_pharmacogenomics_evidence.ipynb)

## Development, versioning, and roadmap

The package is pre-1.0 (`0.4.3` prepared): public APIs and schemas may evolve only under the documented compatibility/deprecation policy. This patch accelerates resumable 00B provisioning without changing the public selector schema or scientific normalization rules. It does not tag, publish, sign, or create a release. Follow the [roadmap](docs/roadmap.md), ADRs, documentation expectations in `AGENTS.md`, and MIT [license](LICENSE).

CI runs these locked gates:

```bash
uv lock --check
uv sync --locked --dev
uv run ruff check .
uv run ruff format --check .
uv run mypy src
GENOME_EVIDENCE_PROFILE=synthetic_ci uv run pytest
GENOME_EVIDENCE_PROFILE=synthetic_ci uv run pytest --nbmake notebooks
uv build
```
