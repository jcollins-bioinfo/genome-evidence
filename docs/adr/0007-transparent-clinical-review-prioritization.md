# ADR-0007: Transparent clinical review prioritization

**Status:** Accepted

## Decision

M4 emits deterministic review bands rather than numeric scores. Routing is controlled by a local, typed, versioned JSON policy whose exact bytes and parsed configuration are checksummed. Bands allocate manual-review attention; they are not pathogenicity, risk, confidence, urgency, penetrance, or actionability.

SCV submissions and VCV aggregates remain independent source assertions. Exact terms and review statuses remain attributed source metadata; review status is not probability. M4 neither votes among submissions nor implements ACMG/AMP criteria.

The initial context is explicitly `germline_constitutional`. Germline assertions can route records; somatic clinical-impact and oncogenicity assertions remain visible but context-only because tumor/sample modeling is absent.

## Consequences

Phenotype fit, family history, inheritance, penetrance, actionability, segregation, and clinical confirmation remain unresolved. Called reference, absent canonical genotype, and discordant calls remain distinct. Absence from the queue is neither comprehensive assay coverage nor a negative result. This conservatism costs convenience but preserves provenance and prevents source vocabulary from becoming a project clinical conclusion.
