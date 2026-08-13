# ADR-0008: Reference-panel population-structure projection

**Status:** Accepted

## Decision

M5 implements fixed-reference PCA projection rather than ancestry percentages. It reports coordinates, source-attributed subset and opaque-sample distances, empirical support envelopes, and leave-one-chromosome-out marker-set sensitivity—never membership, identity, mixture proportions, or clinical conclusions.

Only exact M2 GRCh38 autosomal biallelic-SNV allele keys are accepted. Missing markers are not reference calls; discordant duplicates are excluded without voting. Missing-marker-aware least squares is transparent and avoids implicit mean imputation, but assumes the supplied fixed axes and standardized-genotype model are appropriate. Explicit overlap, chromosome, loading-energy, rank, conditioning, and finite-number gates prevent fabricated coordinates.

Axis signs/scales are bundle-specific. Reference estimation, ascertainment and projection shrinkage—especially in later PCs—remain limitations. Leave-one-chromosome-out variation measures marker-set sensitivity, not statistical uncertainty or assay validation. Reference labels remain exact external descriptors because sampled subsets are not biological essences.

Production bundle construction and real panel data are deferred for licensing and separate scientific review. M5 is non-clinical and non-causal. See [the implementation guide](../ancestry/population-structure.md) for authoritative sources.
