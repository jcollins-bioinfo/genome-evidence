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

Long-running operational state is separate from completed scientific resources:

- `logs/notebooks/00b/<run-key>/events.jsonl` is the append-only, structured execution
  log mirrored to notebook stdout;
- `cache/downloads/normalization/v1/<run-key>/checkpoint.json` binds the checkpoint to
  the selected source and pinned resource identities;
- `cache/downloads/normalization/v1/<run-key>/dbsnp/<assembly>/<query-key>/` contains
  checksum-verified aggregate and per-batch query checkpoints; and
- `cache/downloads/normalization/v1/<run-key>/fasta/ucsc-hg38-gca_000001405.15/`
  retains `hg38.fa.gz` or its resumable `.part` file;
- `references/manifests/normalization/<bundle>.COMPLETED.json` commits a verified
  immutable resource bundle only after every listed artifact and provenance file is
  durable; and
- `config/normalization_resources.PUBLISHING.json` identifies an interrupted selector
  copy, while `config/normalization_resources.COMPLETED.json` is the last-written,
  checksum-bound selector completion marker.

The 24-character run key is deterministically derived from the private source checksum
and pinned dbSNP, Kent, and FASTA identities. Checkpoint directories are therefore not
reused across incompatible sources or pinned-resource changes. They remain private even
though the retrieved dbSNP rows are public. Source-derived identifier-list files exist
only in ephemeral Colab storage while a batch runs; Drive receives verified extracts and
aggregate completion metadata, not the identifier lists.

The GRCh38 archive is fetched from UCSC's fixed `hg38.fa.gz` endpoint and must match its published MD5 `1c9dcaddfa41027f17cd8f7a82c7293b` before decompression. The default checkpoint retains the roughly 938 MiB compressed archive in addition to approximately 3 GiB for the installed FASTA plus FAI. Plan for roughly 4 GiB of Drive usage plus smaller dbSNP checkpoints and logs; a preserved invalid archive can temporarily require another approximately 938 MiB. The checkpoint cache is evictable after a successful, checksum-valid selection, but removing it forfeits download and query resumption if provisioning must later be repeated.

The pinned Kent v479 `bigBedNamedItems` utility first queries local copies of the pinned hg19/hg38 `dbSnp155Common.bb` files. Durable, identity-bound range segments live below the private normalization checkpoint; a verified assembled copy is materialized in ephemeral `/content` for fast queries. The complete 65/68 GiB files are **never downloaded**. Only common-missing rsIDs are sent to their pinned remote indexes. Common absence is not dbSNP absence: every unresolved identifier requires a successfully completed full-index fallback. Genotype tokens are never written to query lists or transmitted.

Plan for roughly 3.6 GiB of additional Drive checkpoint space for both common files on a GRCh37 run, plus approximately 1.8 GiB per required local common copy while active. Before transfer, 00B reports Drive and local free space. Segment files are the durable copy; local assembled/query files are ephemeral. A temporary assembled durable file can coexist during validation, so retain headroom. Drive rename atomicity is not used as a commit boundary: each segment and final resource has a last-written, checksum-bound completion manifest.

### Progress telemetry

Every major step emits timestamped events to notebook stdout and the JSONL log: source
validation, checkpoint discovery, Kent acquisition and validation, dbSNP batch planning,
attempts, retries and adaptive splits, extract parsing, archive download and checksum
verification, decompression, FASTA indexing, artifact publication, provenance, and final
selection publication. Each JSONL line contains a UTC timestamp, severity, stable event
name, human-readable message, and applicable structured fields. The log is mode `0600`
where the filesystem honors POSIX permissions.

Rates have stage-specific meanings:

- direct downloads, hashing, decompression, indexing, and copies report bytes processed,
  human-readable byte rates, percentage, elapsed time, and ETA when total size is known;
- dbSNP BigBed work reports completed batches, identifiers processed, identifiers per
  second, and ETA; and
- no byte-download rate is claimed for BigBed queries because `bigBedNamedItems`
  performs opaque sparse range reads and does not expose reliable transfer-byte
  telemetry.

The log deliberately excludes genotype calls, raw source rows, individual marker
identifiers, and unsanitized identifier-bearing subprocess diagnostics. It is written to
the private Drive workspace, not uploaded to an external logging service.

### Retry and rerun semantics

dbSNP identifiers are deterministically sorted, deduplicated, and queried in bounded
batches. Each successful batch is parsed, checked against its requested identifier set,
hashed, and completed with a manifest. Failed attempts use bounded backoff against the
pinned UCSC endpoint and may be bisected into smaller batches.
Partial output from a failed process is never classified as a completed query. An
identifier absent from a successfully completed query may remain unresolved; an
uncompleted query is an operational failure, not biological missingness.

HTTP downloads retain `.part` bytes and request the remaining range on retry. A resumed
response is appended only when the server returns a matching partial-content range;
otherwise that file restarts from byte zero. Completed FASTA archive bytes are verified
against UCSC's published MD5 and a recorded SHA-256 before use. Gzip decompression cannot
safely resume from arbitrary decoder state, so only that transform restarts while the
verified archive remains cached. Completed decompressed FASTA and dbSNP checkpoints are
checksum-verified before reuse.

Bounded network or workspace-I/O exhaustion becomes `ProvisioningIncomplete`; notebook
00B prints `status: incomplete_resumable`, the durable log and checkpoint paths, and the
instruction to rerun the cell. A rerun reconstructs the same run key, validates existing
components, and continues at the first incomplete unit. Configuration conflicts,
unexpected data, checksum mismatches, unsafe paths, or incompatible completed resources
remain errors. Tolerance never means accepting corrupt or scientifically incomplete
resources. Before a bundle completion marker exists, torn workflow-owned output files
may be repaired from their verified checkpoints. After that marker exists, bundle files
are immutable and any mismatch fails closed. Selector staging is similarly repairable
only while its `PUBLISHING` marker exists; the checksum-validated `COMPLETED` marker is
written last. No workflow step assumes that a Drive FUSE rename is an atomic commit.

### Optional 00B controls

Defaults are intended for ordinary Colab execution. Increase batch size only when the
connection and runtime are stable; smaller batches provide finer checkpoints at greater
process overhead.

| Environment variable | Default | Accepted range | Effect |
|---|---:|---:|---|
| `GENOME_EVIDENCE_DBSNP_BATCH_SIZE` | `5000` | 250–25,000 | Maximum rsIDs in an initial Kent query batch |
| `GENOME_EVIDENCE_QUERY_ATTEMPTS` | `4` | 1–10 | Attempts per dbSNP batch before adaptive splitting or resumable pause |
| `GENOME_EVIDENCE_QUERY_TIMEOUT_SECONDS` | `900` | 60–3,600 | Timeout for one Kent query attempt |
| `GENOME_EVIDENCE_DBSNP_WORKERS` | `6` | 1–12 | Global cap for full-index Kent subprocesses; `1` is serial |
| `GENOME_EVIDENCE_COMMON_DOWNLOAD_SEGMENTS` | `8` | 1–12 | Concurrent HTTP ranges across common downloads |

Each persistent fallback worker reuses only its own local `TMPDIR/udcCache`; caches are never shared. Retries remain bounded and jittered, and adaptive splitting preserves completed child checkpoints. Existing v1 batch manifests are validated against their input digest, assembly, URL, tool identity, output hash/size/count, and returned-ID subset before reuse. Invalid checkpoints are ignored without logging identifiers.

The expected speedup is substantial but not guaranteed: local common hits avoid latency-bound random remote reads, while the remaining independent full queries overlap up to the configured bound. These factors are multiplicative in the ordinary performance sense, not mathematically exponential. Actual improvement depends on the common-hit fraction, UCSC/network behavior, and Colab storage throughput; telemetry reports measured rates and never invents byte rates for Kent sparse queries.

The existing selectors still apply: `GENOME_EVIDENCE_SOURCE_SHA256` chooses one imported
source when more than one exists, and `GENOME_EVIDENCE_SOURCE_BUILD` supplies `GRCh37` or
`GRCh38` only after independent verification when the vendor header is insufficient.

Definitions require exact rsID, chromosome, one-based position, and source assembly agreement. 23andMe documents its GRCh37/GRCh38 raw genotypes on the reference plus strand, so accepted source rows receive `orientation=none` with the vendor assertion URL retained in provenance. Non-rsID markers, absent placements, coordinate disagreements, non-SNV alleles, strand/allele incompatibilities, and ambiguous placements remain unresolved. The GRCh37→GRCh38 map uses only same-rsID placements with identical plus-strand REF and compatible ALT; it is deliberately more conservative than a generic interval liftover.

UCSC does not publish a checksum alongside the Kent binary or remote dbSNP BigBed objects. The provisioner therefore pins the binary's versioned release URL, verifies the executable's documented `bigBedNamedItems` usage signature before use, and records the exact executable SHA-256 and extracted BED SHA-256 values. The Kent utility does not report the suite release number at runtime, so the `kent_version` provenance field denotes the pinned URL release rather than a value asserted by the executable. This provides reproducibility of the local result but is not equivalent to an upstream-signed artifact. The FASTA has the stronger published-checksum verification described above.

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
