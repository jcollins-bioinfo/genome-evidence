# ADR-0005: Explicit resources and derived normalization identity

**Status:** Accepted

## Decision

M2 consumes checksum-validated M1 Parquet and uses only explicit, local, versioned marker, reference, and liftover resources. Source observations remain unchanged. Observation references include M1 identity and source line; canonical variants use assembly/chromosome/position/REF/ALT. Ambiguity and failures are first-class mapping outcomes.

## Consequences

Coordinates and observed tokens alone never establish alleles or builds. This requires curated resources and may leave rows unresolved, but preserves provenance and prevents false biological certainty.
