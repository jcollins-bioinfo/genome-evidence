# Changelog

## 0.6.0 - Unreleased

- Add the M9 declared-pedigree schema, structural validation, checksum-bound completed-M2 subject inputs, exact canonical family alignment, explicit Mendelian transmission enumeration, private atomic artifacts, aggregate-only CLI, notebook 08, and scientific/privacy documentation.
- Keep missing, conflicting, unsupported, and inferred evidence out of ordinary autosomal diploid compatibility; relationship assertions are not relatedness verification, and site transmission is not long-range phase.

## 0.5.0 - Unreleased

- Make `bounded_local_v1` the default normalization provisioning policy, add local ClinVar supplementation, explicit unresolved artifacts, an 80% source-placement guard, and policy-bound provenance/selectors.


- Bind dbSNP query checkpoint v2 manifests to explicit, stable resource identities rather
  than runtime-local Kent execution paths. Narrowly validate and migrate compatible v1
  common-query checkpoints from workflow-owned Colab temporary paths.
- Persist content-bound common-indeterminate leaf dispositions, initialize aggregate
  progress from verified work, and propagate every top-level common-batch completion to
  the overall notebook dashboard without counting resumed identifiers as session throughput.
- Treat checksum-valid completed common BigBeds as immutable local objects by default;
  mutable origin headers are still enforced for incomplete segmented transfers.
- Render the in-place Jupyter dashboard as escaped preformatted HTML instead of a quoted
  Python string, without changing the public normalization selection schema.

## 0.4.4 - Unreleased

- Restore top-level common-query aggregate progress and add a monotonic, versioned
  12-stage workflow-completion dashboard with Unicode, ASCII, and plain renderers.
- Separate current-session transfer bytes from retained/reused segments, explicitly
  label assembly and verification, and reuse checksum-valid completed destinations.

## 0.4.3

- Localize deterministic common-dbSNP validation failures to their minimum query leaves, preserve every independently validated common checkpoint, and authoritatively fall back only for common-missing or common-indeterminate identifiers.
- Add privacy-safe validation categories and explicit, internally consistent common, legacy-full, and full-fallback provenance counts without changing the public selection schema.
- Retain compatible v1 query checkpoints and strict, resumable full-index failure semantics.

## 0.4.2

- Replace serial complete-dbSNP extraction with local common-dbSNP-first querying and a
  globally bounded, isolated-cache full-index fallback.
- Add identity-bound segmented HTTP downloads, durable segment completion manifests,
  safe single-stream fallback, concurrent privacy-safe reporting, and deterministic merge.
- Preserve the public selector schema and conservative normalization semantics while
  recording the new algorithm, common identities, outputs, and worker configuration.

This project follows [Keep a Changelog](https://keepachangelog.com/) structure.

## [Unreleased]

### Added

- Notebook 00B and public workspace provisioning APIs for exact-source dbSNP marker
  extraction, checksummed GRCh38 FASTA/FAI installation, conservative cross-build
  maps, canonical Drive placement, durable selection, and provenance manifests.
- Privacy-safe notebook stdout and structured JSONL telemetry for every long-running
  00B phase, including rates, ETA, retries, checkpoint reuse, and final status.
- Source-bound, checksum-validated 00B checkpoints for resumable HTTP downloads,
  bounded dbSNP query batches, adaptive batch splitting, decompressed FASTA reuse, and
  idempotent final publication.
- Last-written bundle and selector completion markers that make torn Drive/FUSE writes
  repairable before commit and immutable, checksum-validated artifacts afterward.
- An operational personal M1→M2 notebook workflow with content-addressed source
  selection, ephemeral computation, verified Drive publication, and compatible-run
  resolution for downstream notebooks.
- Generic checksum-verified M1/M2/M5 workspace publication and completion registries.
- Indexed FASTA access and source-build marker-definition validation for practical,
  fail-closed personal normalization.
- M8 strict local pharmacogenomics bundle validation, exact-observation candidate
  matching, immutable evidence artifacts, private workspace support, CLI, notebook,
  and documentation.
- A deliberate pre-1.0 version and release policy.

### Changed

- Prepare patch version 0.4.1; installed distribution metadata remains the single
  runtime version authority.
- Notebook 01 resolves the source-compatible resource selection persisted by notebook
  00B, while explicit environment selectors remain supported.
- Colab bootstrap rejects mixed in-memory package revisions before switching a checkout.
- Notebook 04 now resolves a registered compatible M2 run automatically and publishes
  M5 output when a reviewed local population bundle is installed.

### Fixed

- Replace the single long-running remote dbSNP query with retryable, resumable batches
  so a transient UCSC range-query failure no longer discards earlier verified work.
- Preserve partial FASTA download bytes across retries, retain the verified archive
  across later-stage failures, surface sanitized subprocess diagnostics, and report an
  explicit resumable-incomplete outcome without weakening checksum or provenance gates.
- Retry explicit incomplete/chunked HTTP reads and reject incompatible partial-content
  responses before they can be promoted into a completed cache file.

## 0.1.0 — historical development baseline

The repository used this version during M0–M7 development. No Git tag, GitHub
release, or package publication is asserted.
