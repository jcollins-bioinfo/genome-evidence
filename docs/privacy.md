# Privacy and safety

Raw personal genomic, phenotype, laboratory, family-history, and medical data must never be committed. Tests must use only clearly labelled synthetic fixtures. Logs should report identifiers, counts, checksums, and error context without dumping genotypes or medical content.

Future private inputs and local databases should preferably live outside the repository. Derived artifacts and reports can reveal personal information even when raw files are absent; treat them as private and store them accordingly. Secrets belong in an untracked environment or secret manager, never source control.

The `.gitignore` is defense in depth, not a substitute for reviewing staged files and history. These practices do **not** establish or claim regulatory compliance.

M1 never copies the raw input into its output directory and normal CLI output omits the input path and genotype records. Explicit output paths outside the repository are preferred. If an output is inside the worktree, ingestion requires Git to report it ignored; `/runs/private/` is the standard ignored pattern. Manifests use a hash-derived logical source identifier by default rather than the filename or full path. Comment metadata is preserved in private `source_metadata.json` and must itself be treated as identifying.

M4 profiles, candidate tables, rationales, exclusions, and reports are private derived genomic artifacts. The CLI applies the same outside-worktree-or-Git-ignored output guard and prints aggregate counts only. Candidate terms, genotypes, conditions, and accessions must not be emitted to normal logs.
