# Privacy and safety

Raw personal genomic, phenotype, laboratory, family-history, and medical data must never be committed. Tests must use only clearly labelled synthetic fixtures. Logs should report identifiers, counts, checksums, and error context without dumping genotypes or medical content.

Future private inputs and local databases should preferably live outside the repository. Derived artifacts and reports can reveal personal information even when raw files are absent; treat them as private and store them accordingly. Secrets belong in an untracked environment or secret manager, never source control.

The `.gitignore` is defense in depth, not a substitute for reviewing staged files and history. These practices do **not** establish or claim regulatory compliance.

M1 never copies the raw input into its output directory and normal CLI output omits the input path and genotype records. Explicit output paths outside the repository are preferred. If an output is inside the worktree, ingestion requires Git to report it ignored; `/runs/private/` is the standard ignored pattern. Manifests use a hash-derived logical source identifier by default rather than the filename or full path. Comment metadata is preserved in private `source_metadata.json` and must itself be treated as identifying.

M4 profiles, candidate tables, rationales, exclusions, and reports are private derived genomic artifacts. The CLI applies the same outside-worktree-or-Git-ignored output guard and prints aggregate counts only. Candidate terms, genotypes, conditions, and accessions must not be emitted to normal logs.

M5 alignments, coordinates, distances, neighbors, support evaluations, sensitivity records, and reports are private derived genomic artifacts. CLI output is aggregate-only; reference validation and projection never download data.

Private workspaces must remain outside Git. Personal sources, reference genotypes, native VCFs, engine logs, caches, JARs, and executed notebook outputs are not repository material. M6 analysis is offline and aggregate-only; public acquisition cannot inspect target directories. Colab is a Google-managed VM and is not appropriate when cloud processing is unacceptable.
## M8 pharmacogenomics

M8 target analysis is offline. Keep checked bundles under
`references/pharmacogenomics`, downloads in the source-specific caches, and private
outputs under `runs/m8_pharmacogenomics`. CLI/notebook output is aggregate-first;
candidate rows remain private. Never transmit genotype rows, manifests, phenotype
evidence, medication information, or private paths to remote services.

## M9 family descriptors and outputs

Pedigree descriptors and M9 artifacts are private even when identifiers are pseudonymous:
family structure and linked genotype evidence can be identifying. Keep descriptors under
private workspace `inputs/families/` and outputs under `runs/m9_family_analysis/`; both
repository-relative paths are ignored. Files are published with restrictive permissions
where supported. CLI and routine notebook output contain aggregates, not subjects,
relationships, variants, or genotypes. Cloud notebooks process data on a cloud-managed
VM; use local offline execution when that exposure is unacceptable. This is not a claim
of HIPAA, GDPR, or other regulatory compliance.
