# Reference-panel population-structure projection

M5 consumes one validated M2 run and one local, versioned, checksummed PCA bundle. It accepts only exact GRCh38 autosomal biallelic-SNV allele keys. Eligible diploid calls become ALT dosage; missing markers stay missing and discordant duplicate calls are excluded without voting.

The fixed marker-by-component loading matrix is orthonormal. For observed markers, standardized dosage is `x_i=(dosage_i-mean_i)/scale_i`, and coordinates solve `argmin_z ||L_obs z-x_obs||²` with the bundle's explicit `rcond`. Marker-count/fraction, chromosome, loading-energy, rank, condition-number, and finite-value gates fail closed. Distances use reference-PC-standardized Euclidean geometry. Empirical per-subset distance quantiles define support envelopes. Leave-one-chromosome-out reruns are marker-set sensitivity, not confidence intervals.

Coordinates, distances and support are model outputs, not population assignments. Labels describe source-attributed sampled subsets, not discrete natural kinds. PCA can reflect LD, relatedness, imbalance, artifacts and outliers; projected later PCs can shrink. Axes are comparable only within the exact bundle. No production reference bundle is shipped.

Sources: [Patterson et al.](https://doi.org/10.1371/journal.pgen.0020190), [PLINK 2 PCA](https://www.cog-genomics.org/plink/2.0/strat), [1000 Genomes](https://doi.org/10.1038/nature15393), [Privé et al.](https://doi.org/10.1093/bioinformatics/btaa520), [National Academies](https://doi.org/10.17226/26902), and [Royal et al.](https://doi.org/10.1016/j.ajhg.2010.03.011).
