# Privacy and safety

Raw personal genomic, phenotype, laboratory, family-history, and medical data must never be committed. Tests must use only clearly labelled synthetic fixtures. Logs should report identifiers, counts, checksums, and error context without dumping genotypes or medical content.

Future private inputs and local databases should preferably live outside the repository. Derived artifacts and reports can reveal personal information even when raw files are absent; treat them as private and store them accordingly. Secrets belong in an untracked environment or secret manager, never source control.

The `.gitignore` is defense in depth, not a substitute for reviewing staged files and history. These practices do **not** establish or claim regulatory compliance.
