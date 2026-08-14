# Changelog

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
