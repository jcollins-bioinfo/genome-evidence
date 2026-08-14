# Pharmacogenomics model and evidence

M8 is a provenance-first research foundation, not a clinical assay or prescribing
engine. A local immutable bundle separates sources, genes/capabilities, exact loci,
named alleles and haplotype-definition evidence, explicit locus constraints,
source-attributed function assertions, diplotype/phenotype rules, and gene–drug
guideline evidence. Each source requires exact version, URL, retrieval time, byte
hash, license/terms, and content fingerprint. Production acquisition is not shipped;
official resources must be reviewed and prepared privately without target data.

The target path validates M2 lineage and hashes, aligns exact `(assembly,
chromosome, position, REF, ALT)` identities, preserves observation references, and
enumerates all compatible unordered allele pairs. An explicit observed REF/REF can
support a source-defined reference allele; absence, no-call, unresolved mapping, or
unassayed positions cannot. Unmodeled variation and duplicate conflicts prevent
resolution. Limits fail closed without truncation.

Phenotype terms appear only where an exact source-versioned rule applies; diplotype
ambiguity remains visible even when terms converge. Guideline links preserve source
vocabulary but generate no medication choice, dose, alert, or recommendation.

References: [PharmVar](https://www.pharmvar.org/),
[ClinPGx](https://www.clinpgx.org/), [CPIC data releases](https://github.com/cpicpgx/cpic-data/releases),
and [PharmCAT methods and caveats](https://pharmcat.clinpgx.org/Disclaimers/).
