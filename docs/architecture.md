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

Normalization validates M1 manifests and checksums, then uses local checksummed marker definitions, reference sequence, and optional liftover mappings. It writes independent Parquet entities for variants, mappings, genotypes, and candidates; these artifacts are reconstructible analytical inputs, not claimed DuckDB persistence.
