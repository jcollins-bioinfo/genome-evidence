# ADR-0001: Immutable observations and versioned derivations

**Status:** Accepted

## Decision

Preserve source observations immutably. Normalization, liftover, correction, annotation, and interpretation create separately identified, versioned records with run and input provenance. DuckDB provides relational query state that can be reconstructed from source and derived artifacts where practical.

## Rationale and consequences

Destructive annotation loses what was measured, makes changing references difficult to audit, and confuses source with interpretation. Immutability costs storage and requires explicit joins, but enables reproducibility, correction, comparison, and reprocessing.
