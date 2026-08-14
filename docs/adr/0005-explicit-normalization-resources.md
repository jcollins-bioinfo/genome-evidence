# ADR-0005: Explicit resources and derived normalization identity

**Status:** Accepted

## Decision

M2 consumes checksum-validated M1 Parquet and uses only explicit, local, versioned marker, reference, and liftover resources. M0B may acquire those resources separately from M2, but must bind them to one source checksum, retain upstream URLs/releases/checksums, and persist only exact source-build rsID/coordinate matches. Marker definitions must match the resolved source assembly. Large FASTA references require a checksummed adjacent FAI for bounded random access. Source observations remain unchanged. Observation references include M1 identity and source line; canonical variants use assembly/chromosome/position/REF/ALT. Ambiguity and failures are first-class mapping outcomes.

## Consequences

Coordinates and observed tokens alone never establish alleles or builds. This requires curated resources and may leave rows unresolved, but preserves provenance and prevents false biological certainty.
