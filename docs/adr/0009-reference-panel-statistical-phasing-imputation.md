# ADR 0009: local reference-panel statistical phasing and imputation

- Status: accepted infrastructure; execution incomplete pending mandatory real-backend smoke

Use pinned external Beagle 5.5 rather than a home-grown HMM. Analyze only GRCh38 autosomes and exact canonical alleles from M2; M5 neither selects nor subsets the source-agnostic reference panel. Preserve observations and represent statistical phase and imputation as separate inferences. Haplotype labels are arbitrary: do not claim parent of origin or invent phase blocks, PS, GP, or confidence.

Runs are per chromosome and identities bind upstream, code, configuration, tool, panel and target hashes. Drive completion copies to a unique path, re-hashes, and writes `COMPLETED.json` last. `personal_drive` fails closed on missing private prerequisites; `synthetic_ci` is fabricated, deterministic and offline. Masking measures only internal consistency and cannot establish external or clinical validity. M6 has no automatic downstream interface.
