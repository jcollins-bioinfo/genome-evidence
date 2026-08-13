# Epistemic contract

This contract governs domain semantics. Convenience must never erase how a claim became known.

## Terms

- **Observation:** an immutable record of what an assay or person supplied, including original representation, method, sample, and call status. It is not silently corrected.
- **Inference:** a statistically or logically derived state, identified with its method, inputs, quality, uncertainty, and run provenance. It is not a measurement.
- **External assertion:** a claim made by a named, versioned outside source, with accession and retrieval date; it is not timeless truth.
- **Evidence:** a typed, attributable item that may support, oppose, or contextualize a proposition. Heterogeneous evidence does not become interchangeable merely by aggregation.
- **Model output:** the result of an identified model and configuration. It is conditional on assumptions and does not automatically constitute a diagnosis or probability applicable to this subject.
- **Missingness:** absence due to no call, lack of assay coverage, filtering, or unavailable data. The reason should be retained where known. Missing is a state, not an allele.
- **Uncertainty:** limits on measurement, mapping, inference, interpretation, or transportability. Assay confidence and interpretation confidence are separate dimensions.
- **Provenance:** the chain connecting a derived record to input records and checksums, software version and git commit, configuration hash, reference versions, transformation/model, and pipeline run.
- **Interpretation status:** a project vocabulary (`established`, `supported`, `provisional`, `speculative`, `indeterminate`) for communicating maturity. It is not by itself a clinically validated grading system.

## Hard invariants

1. Observations, canonical representations, inferences, assertions, phenotypes, and outputs retain distinct types and identifiers.
2. Source observations are immutable. Normalization, strand correction, liftover, annotation, and interpretation create derived records linked by provenance.
3. Every derived result requires run and input provenance. External evidence also records source version and retrieval time.
4. Quantitative clinical precision requires a validated, transportable model. Otherwise report structured evidence, uncertainty, hypotheses, and unresolved questions—not invented probabilities, diagnoses, classifications, or recommendations.

## Invalid transformations

```text
missing locus -> homozygous reference
imputed pathogenic allele -> confirmed pathogenic genotype
high assay confidence -> high interpretation confidence
current external assertion -> timeless biological truth
normalized allele -> overwrite source observation
heterogeneous evidence count -> absolute disease probability
```

A no-call and a called reference/reference genotype must remain queryably different. An imputed allele may motivate confirmation; it cannot be relabeled as directly observed.

## Source-ingestion boundary (M1)

A parsed vendor genotype remains a source token, not a canonical allele or biological genotype. Lexical categories describe token shape only; they do not infer strand, REF/ALT, ploidy, or biological zygosity. User-provided assembly overrides remain distinguishable from vendor declarations, and coordinates never imply a build.
