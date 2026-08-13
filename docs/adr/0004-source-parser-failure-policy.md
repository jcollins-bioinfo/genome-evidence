# ADR-0004: Strict and lenient source-parser failure policy

**Status:** Accepted

## Decision

23andMe ingestion defaults to strict parsing: a structurally invalid record aborts before artifacts are completed. An explicit lenient mode retains safe line-numbered malformed-record diagnostics and continues. Neither mode repairs a record. All valid records, including duplicates, retain source order. In-worktree output is allowed only at Git-ignored paths.

## Rationale and consequences

Formal analysis should not silently proceed past structural corruption, while forensic inspection needs complete error accounting. Separate modes make that choice explicit. Rejecting unignored in-repository output reduces accidental commits of identifying derived data; callers may always use explicit private paths outside the repository.
