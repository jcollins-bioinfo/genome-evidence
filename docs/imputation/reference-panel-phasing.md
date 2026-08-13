# Local reference-panel phasing and imputation

M6 is autosomal GRCh38 infrastructure for diploid biallelic A/C/G/T SNVs. It depends scientifically only on a completely validated M2 run and joins the local phased panel solely by exact `(assembly, chromosome, position, REF, ALT)`. It never performs rsID/position matching, swapping, complementing, implicit liftover, panel selection from M5, genotype refinement, or network analysis.

Beagle 5.5 `27Feb25.75f` is pinned in `config/beagle-5.5.json`. The JAR and production panel are not distributed. Acquisition is a separate public-resource action; analysis refuses downloads. Observed M2 allele multisets must remain unchanged. Statistical phase (with arbitrary haplotype labels) and imputed genotypes are separate source-attributed inferences; DS is ALT dosage, AP1/AP2 are ALT probabilities for each arbitrary haplotype, and DR2 is estimated squared correlation, not correctness probability. Fields are retained only when declared by the engine header; GP, PS, and confidence are never invented.

Masked-marker results are internal consistency assessment, not external or clinical validation. Results depend on switch error, panel, marker density, MAF and LD; consumer arrays are limited. Absence is not 0/0. Outputs support no clinical, ancestry, identity, or family conclusion and have no automatic M3/M4/PRS/downstream handoff.

The contracts and validators are implemented, but end-to-end publication remains explicitly incomplete until a real Beagle smoke with a sufficiently valid fabricated phased panel yields both preserved observed phase and an imputed marker. No production reference resource ships in this repository.
