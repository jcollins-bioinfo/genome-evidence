# Architecture

Genome Evidence is a provenance-first pipeline with typed boundaries. Adapters write immutable source observations; transformations emit new, versioned artifacts; DuckDB provides reconstructible analytical query state over those artifacts. Pydantic models define API boundaries, Polars supports dataframe transformations, and Parquet is the primary derived/intermediate tabular format.

## Planned layers

1. **Observation** — source-faithful assay, phenotype, laboratory, and family inputs.
2. **Normalization** — canonical alleles, builds, strand handling, and explicit mappings.
3. **Annotation** — versioned external biological assertions.
4. **Phasing/imputation** — probabilistic, quality-bearing inference records.
5. **Evidence integration** — typed relationships without collapsing epistemic categories.
6. **Mendelian analysis** — inheritance-aware reasoning.
7. **Polygenic scoring** — explicit model/version and portability limits.
8. **Pharmacogenomics** — modelled diplotypes and evidence, not automatic advice.
9. **Ancestry** — scoped statistical estimates.
10. **Family analysis** — relationships and multiple subjects/samples.
11. **Phenotype integration** — measured, reported, and family evidence remain distinct.
12. **Hypothesis/risk models** — assumptions, uncertainty, and transportability are explicit.
13. **Value-of-information** — identify useful confirmation or missing inputs.
14. **Reporting** — evidence-oriented views without synthetic clinical precision.

M1 ingestion and M2 resource-driven normalization are implemented. A future WGS/VCF/BCF adapter enters through the immutable observation layer, so downstream canonical mappings and evidence models do not need redesign. 23andMe is one source adapter, not the conceptual model.

## Confidence and provenance

Measurement/mapping quality, statistical inference quality, and interpretation support are independent. Derived records reference input IDs and a pipeline run containing input hashes, configuration hash, software/git identity, transformation, and reference versions. Raw records are never destructively enriched.

## M1 implementation

The observation layer now includes a source-faithful 23andMe adapter with strict/lenient structural parsing and a separate descriptive QC package. The adapter processes logical lines from one exact-file byte read, retains source order/line numbers and vendor tokens, and emits provenance-linked JSON, Parquet, and Markdown artifacts. It creates no `Variant` or `ObservationVariantMapping`; those remain M2 responsibilities.

## M2 implementation

Normalization validates M1 manifests and checksums, then uses source-build-matched marker definitions, a checksummed local reference sequence (FAI-backed for large FASTA), and optional liftover mappings. M0B provisions these independently of genotype calls. Its default `bounded_local_v1` plan downloads identity-bound, resumable segments of each pinned Common and ClinVar BigBed, queries assembled local copies, and queries ClinVar only for Common-unresolved identifiers. It never contacts the complete remote indexes. Local-track non-return is explicit unresolved coverage—not scientific absence—and is published when the 80% gross exact-placement guard passes. Deterministic merge and the existing exact placement/SNV validation preserve scientific semantics regardless of completion order or worker count. The FASTA is checksum-verified, and selector-last publication remains fail-closed.

## M3 implementation

The offline ClinVar adapter streams a checksummed VCV XML snapshot into typed source snapshot, variant representation, assertion, relationship, and condition records. Linking validates every upstream artifact and compares only the exact canonical allele tuple. Assertions and links remain separate from observations, genotypes, and interpretations.

## M4 implementation

Prioritization consumes three mutually consistent, checksum-validated upstream manifests and a local typed policy. It builds a profile for every exact-linked M2 variant, computes genotype evidence without voting discordant rows, and writes atomic Parquet/JSON/Markdown artifacts. Structured candidate-to-assertion and rationale tables preserve every routing input. Runtime identity is distinct from deterministic scientific record identity.

## M5 implementation

Population structure consumes only checksum-valid M2 artifacts and a local strict PCA reference bundle. Exact allele alignment, partial-marker least squares, numerical gates, model-dependent distances/support, and marker-set sensitivity remain a separate non-clinical model-output layer.

## Private workspace and M6 boundary

The `workspace` package owns provider-neutral durable paths, content-addressed source import, relative JSON run registries, and M1/M2/M5 publication. Registered resolution validates the completion marker and every declared artifact before returning a path. M6 consumes validated M2 directly through a local-only engine adapter and source-agnostic checked reference contract. Its phase and imputation tables are inference stores, never observation stores or automatic downstream inputs. Publication uses copy, destination re-hashing, then a last-written completion marker rather than filesystem rename atomicity.

## M8 pharmacogenomics

The pharmacogenomics package validates a local immutable evidence bundle, validates
an upstream M2 run, creates exact observed-locus evidence, enumerates bounded
candidate diplotypes, and atomically publishes typed artifacts. Acquisition is
separate from target analysis; unsupported methods fail closed. See ADR 0011.
