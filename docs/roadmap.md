# Roadmap

- **M0 — Bootstrap:** epistemic contract, domain boundaries, provenance, storage scaffold, CLI, and quality gates.
- **M1 — 23andMe ingestion/QC (implemented):** source-faithful adapter and assay-level quality reporting.
- **M2 — Canonical normalization (implemented):** explicit allele mappings, strand handling, and liftover provenance.
- **M3 — External annotation/evidence (implemented):** verified in the locked M4 hardening CI history.
- **M4 — Clinical variant prioritization (implemented):** verified in the locked M4 hardening CI history. Evidence-oriented manual-review routing only, without automatic classification.
- **M5 — Reference-panel population structure (implemented infrastructure):** exact M2 alignment, PCA projection, distances, support, and marker-set sensitivity. Synthetic demonstration only; no production reference bundle is shipped.
- **M6 — Phasing/imputation (incomplete infrastructure):** workspace, exact-input/reference/tool contracts, and offline adapter are implemented; mandatory real-Beagle imputation smoke and completed-result publication remain incomplete. No production resources ship.
- **M7 — Polygenic scoring (implemented foundation):** versioned models and portability limitations.
- **M8 — Pharmacogenomics (implemented foundation):** strict local evidence bundles and conservative star-allele candidates; no production bundle, structural calling, clinical validation, or recommendations.
- **M9 — Family-aware analysis:** relationships, segregation, and inheritance evidence.
- **M10 — Phenotype integration:** measured and reported phenotype inputs.
- **M11 — Hypothesis engine:** structured hypotheses and unresolved questions.
- **M12 — Value-of-information:** prioritize confirmatory measurements and missing evidence.
- **M13 — Reporting/UI:** privacy-aware evidence reporting.
- **M14 — WGS ingestion:** scalable VCF/BCF/WGS observation adapters.

Later milestones are proposals, not implemented functionality or clinical claims. Each requires review against the epistemic contract.
