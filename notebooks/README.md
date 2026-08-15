# Notebook launch index

The locked local environment and `synthetic_ci` nbmake CI are canonical. `personal_drive` uses checked private resources and never falls back to synthetic data. Colab runs on a Google-managed VM; use the local path when cloud processing is unacceptable. Opening a badge does not auto-run, mount Drive, grant access, create folders, or guarantee the lockfile environment.

## Clean-runtime behavior

- Notebook 00 can bootstrap a fresh Colab kernel, initialize the private workspace, and import one explicitly selected 23andMe inbox file.
- Notebook 00B downloads and validates the source-compatible dbSNP marker subset, GRCh38 FASTA/FAI, and any required variant-specific cross-build map into canonical Drive subdirectories. It streams progress to the cell and a private JSONL log, checkpoints downloads and dbSNP batches, and persists provenance and selectors for later kernels only after every required artifact verifies.
- Notebook 01 consumes the content-addressed source plus the durable 00B resource selection. It computes ephemerally, publishes checksum-verified M1/M2 runs to Drive, and registers the latest compatible M2 run.
- Notebooks 02–03 still execute fabricated demonstrations in personal mode; they do not yet consume private production inputs.
- Notebook 04 automatically resolves the latest compatible M2 run and publishes M5 output when exactly one reviewed population-reference bundle is installed.
- Notebooks 05–07 stop when their required completed-run, reference-bundle, model, gene, or tool prerequisites are absent. No personal path substitutes synthetic resources.
- CI executes the notebooks under the preinstalled `synthetic_ci` profile. The separate personal-bootstrap regression test covers source installation and same-process import in a fresh interpreter.

## 00B progress and resumption

00B uses **common-first + bounded-parallel fallback**. It downloads pinned common
BigBeds with 8 resumable segments by default, queries verified local copies, then sends
only common-missing and common-indeterminate identifiers to the authoritative complete remote indexes with 6 isolated-cache Kent
workers. Configure `GENOME_EVIDENCE_COMMON_DOWNLOAD_SEGMENTS` and
`GENOME_EVIDENCE_DBSNP_WORKERS` from 1–12; use `1` for serial behavior. The complete
65/68 GiB BigBeds are never downloaded. Restart the runtime and rerun the same cell after
an interruption—do not delete the checkpoint.

During personal execution, 00B prints timestamped phase transitions, retries, aggregate
counts, rates, and ETA directly below the provisioning cell. The same events are appended
as structured JSON Lines to
`logs/notebooks/00b/<run-key>/events.jsonl` in the private workspace. Logs contain
aggregate marker counts, public resource metadata, operational paths, checksums, and
diagnostics; they do not contain genotype calls, raw source rows, or individual marker
identifiers.

Durable state is stored under
`cache/downloads/normalization/v1/<run-key>/`. A rerun validates each completion marker
and checksum before reusing it. Successful dbSNP batches and adaptive-split children are not repeated. A deterministic common-output rejection is categorized without identifiers, localized to its minimum leaf, and routed to full fallback without discarding validated siblings. An
interrupted HTTP transfer continues from its retained `.part` bytes when the server
honors the requested byte range. A partial gzip decompression restarts that transform
from the retained, checksum-verified archive because gzip decoder state cannot safely be
continued at an arbitrary output offset.

If bounded network retries are exhausted, the cell prints
`status: incomplete_resumable` plus its log and checkpoint paths instead of discarding
verified work. Rerun the same provisioning cell. Configuration or integrity failures
remain fail-closed: 00B does not publish `config/normalization_resources.json` until all
required artifacts and provenance validate.

| # | Role | Prerequisite | Repository | Colab |
|---:|---|---|---|---|
| 00 | initialize private workspace | None | [00_initialize_private_workspace.ipynb](00_initialize_private_workspace.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/00_initialize_private_workspace.ipynb) |
| 00B | provision normalization resources | 00 and one imported source | [00b_provision_normalization_resources.ipynb](00b_provision_normalization_resources.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/00b_provision_normalization_resources.ipynb) |
| 01 | ingest, normalize, publish M1/M2 | 00B completed resource selection | [01_ingest_and_normalize_genome.ipynb](01_ingest_and_normalize_genome.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/01_ingest_and_normalize_genome.ipynb) |
| 02 | versioned external evidence | compatible completed M2 and local evidence | [02_versioned_external_evidence.ipynb](02_versioned_external_evidence.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/02_versioned_external_evidence.ipynb) |
| 03 | evidence-oriented clinical review | compatible completed M3 | [03_evidence_oriented_clinical_review.ipynb](03_evidence_oriented_clinical_review.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/03_evidence_oriented_clinical_review.ipynb) |
| 04 | population structure projection | compatible completed M2 and reviewed M5 bundle | [04_population_structure_projection.ipynb](04_population_structure_projection.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/04_population_structure_projection.ipynb) |
| 05 | reference-panel phasing and imputation | compatible completed M2, pinned engine, and reviewed M6 bundle | [05_reference_panel_phasing_and_imputation.ipynb](05_reference_panel_phasing_and_imputation.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/05_reference_panel_phasing_and_imputation.ipynb) |
| 06 | versioned polygenic-score calculation | compatible completed M2 and explicit checked local PGS bundle | [06_polygenic_score_calculation.ipynb](06_polygenic_score_calculation.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/06_polygenic_score_calculation.ipynb) |
| 07 | pharmacogenomic candidate evidence foundations | compatible completed M2, checked local PGx bundle, and explicit genes | [07_pharmacogenomics_evidence.ipynb](07_pharmacogenomics_evidence.ipynb) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/notebooks/07_pharmacogenomics_evidence.ipynb) |
