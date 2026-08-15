# ADR 0012: bounded local dbSNP coverage is the provisioning default

## Status

Accepted for 0.5.0.

## Context

Sparse name-index queries against the complete UCSC dbSNP 155 BigBeds (about 65 GiB
for hg19 and 68 GiB for hg38) have remote, per-identifier latency. More Kent workers
can increase CDN and UDC-cache contention; concurrency cannot bound that external
latency. This makes the former mandatory fallback unsuitable for a Colab default.

A missing local definition is unresolved coverage, not evidence that a variant is
absent, invariant, reference, or safely inferable. M2 already preserves this boundary
as `MARKER_DEFINITION_ABSENT`.

## Decision

`bounded_local_v1` is selected by default through
`GENOME_EVIDENCE_DBSNP_COVERAGE_POLICY`. For each required assembly, provisioning
checksum-verifies and queries a local `dbSnp155Common.bb`, then queries a local
`dbSnp155ClinVar.bb` only for Common-unresolved canonical rsIDs. It deterministically
merges the tracks and never contacts or downloads the complete BigBeds. Conflicting
placements are unresolved rather than selected by track priority. `dbSnp155Mult.bb`
is not an acceptance source.

Publication is allowed with residual unresolved markers when at least 80% of requested
canonical rsIDs have exact source-build placements. Falling below that gross sanity
guard fails closed as a likely assembly, corruption, or compatibility error. The
private checksummed unresolved artifact records reasons without genotype calls.

`full_remote_v1` remains an explicit advanced legacy opt-in. It can be slow and is
potentially unsuitable for Colab. Its checkpoints do not determine bounded-local
output. The established run key remains unchanged, allowing Common query/download,
FASTA, and FAI checkpoints to be reused across ephemeral execution paths. Old
selectors and immutable bundles remain readable, but only a selector carrying the
requested policy and current builder identity satisfies provisioning.

## Consequences

For a GRCh37 source the default footprint is approximately two 1.8 GiB Common files,
two 75 MiB ClinVar files, the GRCh38 FASTA archive and decompressed FASTA/FAI, plus
headroom for temporary local copies. It requires no 65/68 GiB complete BigBed. Once
caches exist, work is bounded by local sequential reads and deterministic
transformations; no exact runtime is promised. Full enrichment is demand-driven or
performed in a controlled offline environment.
