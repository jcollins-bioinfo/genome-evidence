# ADR-0007: Transparent clinical review prioritization

**Status:** Accepted

## Decision

M4 emits deterministic review bands rather than numeric scores. Routing is controlled by a local, typed, versioned JSON policy whose exact bytes and parsed configuration are checksummed. Bands allocate manual-review attention; they are not pathogenicity, risk, confidence, urgency, penetrance, or actionability.

SCV submissions and VCV aggregates remain independent source assertions. Exact terms and review statuses remain attributed source metadata; review status is not probability. M4 neither votes among submissions nor implements ACMG/AMP criteria.

The initial context is explicitly `germline_constitutional`. Germline assertions can route records; somatic clinical-impact and oncogenicity assertions remain visible but context-only because tumor/sample modeling is absent.

## Consequences

Phenotype fit, family history, inheritance, penetrance, actionability, segregation, and clinical confirmation remain unresolved. Called reference, absent canonical genotype, and discordant calls remain distinct. Absence from the queue is neither comprehensive assay coverage nor a negative result. This conservatism costs convenience but preserves provenance and prevents source vocabulary from becoming a project clinical conclusion.

## Terminology sources

The routing vocabulary is source-attributed and follows the boundaries described by the
[ACMG/AMP interpretation framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/),
[ClinGen variant-classification guidance](https://www.clinicalgenome.org/tools/clingen-variant-classification-guidance/),
[ClinVar classification terminology](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/), and
[ClinVar review-status definitions](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
M4 does not implement ACMG/AMP criteria. Clinical actionability is a separate ClinGen curation
activity, and direct-to-consumer observations retain the limitations described by
[FDA guidance](https://www.fda.gov/medical-devices/in-vitro-diagnostics/direct-consumer-tests).
