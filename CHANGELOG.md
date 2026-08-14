# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/) structure.

## [Unreleased]

### Added

- Notebook 00B and public workspace provisioning APIs for exact-source dbSNP marker
  extraction, checksummed GRCh38 FASTA/FAI installation, conservative cross-build
  maps, canonical Drive placement, durable selection, and provenance manifests.
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

- Prepare package version 0.4.0; installed distribution metadata remains the single
  runtime version authority.
- Notebook 01 resolves the source-compatible resource selection persisted by notebook
  00B, while explicit environment selectors remain supported.
- Colab bootstrap rejects mixed in-memory package revisions before switching a checkout.
- Notebook 04 now resolves a registered compatible M2 run automatically and publishes
  M5 output when a reviewed local population bundle is installed.

## 0.1.0 — historical development baseline

The repository used this version during M0–M7 development. No Git tag, GitHub
release, or package publication is asserted.
