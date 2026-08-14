# ADR 0011: pharmacogenomics evidence and star-allele foundations

## Decision

M8 accepts only checksum-valid M2 observed canonical genotypes and a strict local
`genome-evidence-pgx-bundle/v1`. Exact GRCh38 canonical allele identity is the only
join. The generic O(A²L) matcher retains all compatible unordered diplotypes,
reports missingness and ambiguity, and never ranks candidates. Package, schemas,
milestone, and source releases remain independent.

Nomenclature/haplotype evidence, function assertions, phenotype rules, guideline
links, software candidates, clinical laboratory results, and treatment decisions are
separate layers. Guideline records are evidence links, never recommendation prose.

## Scope and consequences

Only explicitly declared autosomal, diploid, normalized biallelic small-variant star
definitions execute. CYP2D6, CNV/SV/hybrids, HLA, mitochondrial, sex-chromosome,
special-algorithm, and unsupported-ploidy paths fail closed. Missing observations
never imply reference or `*1`; even one surviving pair remains insufficient when a
required locus is missing. No production bundle is redistributed.
