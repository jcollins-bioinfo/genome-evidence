# Canonical private workspace

The provider-neutral workspace defaults in personal Colab to `/content/drive/MyDrive/genome-evidence-private`; `GENOME_EVIDENCE_WORKSPACE` overrides it. Manual uploads belong only in `inputs/inbox/23andme/`. Initialization is idempotent and rejects Git checkouts. Exact bytes are copied to `inputs/raw/23andme/<sha256>/genome.txt`; cache is evictable while inputs, references, registries, runs, reports, and exports are durable.

Use only pseudonyms such as `subject-0001`. The JSON registries contain workspace-relative paths; no symlink is used for latest. Compute in ephemeral local storage, copy into a unique Drive destination, verify every hash there, then write `COMPLETED.json` last. Incomplete runs remain under `_incomplete` and are not registered. Drive FUSE rename atomicity is never assumed.

Colab executes on a Google-managed VM. Users who require no cloud processing must run the same path locally. Personal mode never substitutes synthetic data, and production marker, FASTA, evidence, M5, and M6 bundles must be installed and validated explicitly.

## Notebook 01 normalization resources

Notebook 01 discovers exactly one regular file in each applicable canonical directory. If a directory contains multiple versions, select one with the corresponding environment variable; relative values are resolved from the workspace and all selected resources must remain inside it.

| Resource | Canonical directory | Explicit selector |
|---|---|---|
| source-build marker definitions JSON | `references/markers/23andme/` | `GENOME_EVIDENCE_MARKER_DEFINITIONS` |
| GRCh38 FASTA | `references/genome/grch38/` | `GENOME_EVIDENCE_GRCH38_FASTA` |
| GRCh37→GRCh38 JSON liftover map, when required | `references/liftover/grch37_to_grch38/` | `GENOME_EVIDENCE_GRCH37_TO_GRCH38_LIFTOVER` |

A large FASTA requires its standard adjacent `.fai` index; both FASTA and index identities are recorded. Marker definitions must declare the same assembly as the source. If the source file lacks verified build metadata, set `GENOME_EVIDENCE_SOURCE_BUILD` to `GRCh37` or `GRCh38` only after establishing it independently. When multiple imported sources exist, select one by content ID with `GENOME_EVIDENCE_SOURCE_SHA256`.

Notebook 01 never downloads or fabricates these resources. After successful computation it registers `registry/latest/M2.json` and also sets `GENOME_EVIDENCE_NORMALIZATION_RUN` for the current kernel. Later notebooks can resolve the registered run without relying on a process environment variable.

## Notebook 04 population reference

Install exactly one reviewed bundle below `references/population_structure/<bundle>/` or set `GENOME_EVIDENCE_POPULATION_BUNDLE` to a workspace-contained bundle directory. The repository does not ship a production population reference. With a valid bundle, notebook 04 resolves the latest compatible GRCh38 M2 run, performs the offline projection, publishes a verified M5 run, and registers it.
