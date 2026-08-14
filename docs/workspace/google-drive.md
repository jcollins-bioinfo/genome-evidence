# Canonical private workspace

The provider-neutral workspace defaults in personal Colab to `/content/drive/MyDrive/genome-evidence-private`; `GENOME_EVIDENCE_WORKSPACE` overrides it. Manual uploads belong only in `inputs/inbox/23andme/`. Initialization is idempotent and rejects Git checkouts. Exact bytes are copied to `inputs/raw/23andme/<sha256>/genome.txt`; cache is evictable while inputs, references, registries, runs, reports, and exports are durable.

Use only pseudonyms such as `subject-0001`. The JSON registries contain workspace-relative paths; no symlink is used for latest. Compute in ephemeral local storage, copy into a unique Drive destination, verify every hash there, then write `COMPLETED.json` last. Incomplete runs remain under `_incomplete` and are not registered. Drive FUSE rename atomicity is never assumed.

Colab executes on a Google-managed VM. Users who require no cloud processing must run the same path locally. Personal mode never substitutes synthetic data. Notebook 00B installs the M2 marker, FASTA, and cross-build resources described below; production evidence, M5, M6, M7, and M8 bundles remain separately reviewed inputs.

## Notebook 00B resource provisioning

Run notebook 00B after notebook 00 has imported exactly one content-addressed source. It writes:

- `references/markers/23andme/dbsnp155-<build>-<bundle>/marker-definitions.json` plus the exact source and GRCh38 dbSNP extracts;
- `references/genome/grch38/ucsc-hg38-gca_000001405.15/hg38.fa` and adjacent `hg38.fa.fai`;
- `references/liftover/grch37_to_grch38/dbsnp155-<bundle>/variant-coordinate-map.json` for GRCh37 sources;
- `references/manifests/normalization/<bundle>.json` with upstream identities, retrieval time, artifact hashes/sizes, and aggregate resolution counts; and
- `config/normalization_resources.json`, a source-checksum-bound set of durable workspace-relative selectors consumed automatically by notebook 01.

The GRCh38 archive is fetched from UCSC's fixed `hg38.fa.gz` endpoint and must match its published MD5 `1c9dcaddfa41027f17cd8f7a82c7293b` before decompression. Expect roughly 938 MB of network transfer and approximately 3 GB of durable Drive storage for FASTA plus FAI. The pinned Kent v479 `bigBedNamedItems` utility queries fixed UCSC dbSNP 155 hg19/hg38 BigBed indexes by rsID; genotype tokens are never written to the identifier list or transmitted. HTTPS range requests retrieve only indexed records needed for the source.

Definitions require exact rsID, chromosome, one-based position, and source assembly agreement. 23andMe documents its GRCh37/GRCh38 raw genotypes on the reference plus strand, so accepted source rows receive `orientation=none` with the vendor assertion URL retained in provenance. Non-rsID markers, absent placements, coordinate disagreements, non-SNV alleles, strand/allele incompatibilities, and ambiguous placements remain unresolved. The GRCh37→GRCh38 map uses only same-rsID placements with identical plus-strand REF and compatible ALT; it is deliberately more conservative than a generic interval liftover.

UCSC does not publish a checksum alongside the Kent binary or remote dbSNP BigBed objects. The provisioner therefore pins their release/version URLs, validates the executable-reported Kent version, and records the exact executable SHA-256 and extracted BED SHA-256 values. This provides reproducibility of the local result but is not equivalent to an upstream-signed artifact. The FASTA has the stronger published-checksum verification described above.

## Notebook 01 normalization resources

Notebook 01 first loads `config/normalization_resources.json` when it matches the selected source checksum. Explicit environment variables override the durable selection. Without either, it falls back to discovering exactly one regular file in each applicable canonical directory. Relative values are resolved from the workspace and all selected resources must remain inside it.

| Resource | Canonical directory | Explicit selector |
|---|---|---|
| source-build marker definitions JSON | `references/markers/23andme/` | `GENOME_EVIDENCE_MARKER_DEFINITIONS` |
| GRCh38 FASTA | `references/genome/grch38/` | `GENOME_EVIDENCE_GRCH38_FASTA` |
| GRCh37→GRCh38 JSON liftover map, when required | `references/liftover/grch37_to_grch38/` | `GENOME_EVIDENCE_GRCH37_TO_GRCH38_LIFTOVER` |

A large FASTA requires its standard adjacent `.fai` index; both FASTA and index identities are recorded. Marker definitions must declare the same assembly as the source. If the source file lacks verified build metadata, set `GENOME_EVIDENCE_SOURCE_BUILD` to `GRCh37` or `GRCh38` only after establishing it independently. When multiple imported sources exist, select one by content ID with `GENOME_EVIDENCE_SOURCE_SHA256`.

Notebook 01 never downloads or fabricates resources; acquisition remains isolated in 00B. After successful computation it registers `registry/latest/M2.json` and also sets `GENOME_EVIDENCE_NORMALIZATION_RUN` for the current kernel. Later notebooks can resolve the registered run without relying on a process environment variable.

## Notebook 04 population reference

Install exactly one reviewed bundle below `references/population_structure/<bundle>/` or set `GENOME_EVIDENCE_POPULATION_BUNDLE` to a workspace-contained bundle directory. The repository does not ship a production population reference. With a valid bundle, notebook 04 resolves the latest compatible GRCh38 M2 run, performs the offline projection, publishes a verified M5 run, and registers it.
