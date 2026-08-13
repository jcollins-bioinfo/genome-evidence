# ADR-0006: Snapshot-bound assertions and exact allele links

**Status:** Accepted

## Decision

M3 identifies a ClinVar VCV XML snapshot from authoritative XML release metadata and the SHA-256 of the exact local file. SCV submissions and VCV aggregates are separate immutable assertions connected by explicit relationships. Each assertion has a logical accession/version key, a content fingerprint, and a snapshot-specific instance ID.

External variant representations link to M2 variants only when assembly, canonical chromosome, one-based position, REF, and ALT are exactly equal. The link is separate from every assertion. M3 performs no normalization, conflict resolution, interpretation, or carrier determination.

## Consequences

Conflicts and changed content remain detectable and queryable across snapshots. Some valid ClinVar structures remain unsupported or unmatched rather than being guessed. A VCV annotation can exist for an M2 variant derived from a no-call or reference/reference observation; it never claims that ALT was observed.
