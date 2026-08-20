# ADR 0013: declared pedigrees and site-local transmission enumeration

- Status: accepted for the M9 foundation

M9 treats biological-parent edges as source-attributed user assertions, not genetically
verified relationships. It consumes only checksum-valid completed M2 direct canonical
observations and enumerates allele transmission for autosomal biallelic diploid trios and
duos. Missing, conflicting, unsupported, and inferred evidence cannot be coerced into
calls. Compatibility and transmission uniqueness are separate typed axes.

This conservative design makes an inconsistency an auditable conditional result without
choosing mutation, non-parentage, identity, normalization, or assay-error explanations.
Transmission assignments apply to one site and declared edges only. M9 does not perform
kinship, de novo calling, parent-of-origin inference beyond that context, long-range
phase, recombination, linkage, penetrance, pathogenicity, or clinical interpretation.
