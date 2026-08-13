# ADR-0002: Observation versus inference

**Status:** Accepted

## Decision

Direct observations and genotype inferences are distinct immutable domain types and tables. No convenience conversion may relabel an inference as an observation. Missing/no-call states cannot produce reference calls without new evidence. Mapping records connect observations to canonical alleles without rewriting either.

## Rationale and consequences

The separation prevents probabilistic or transformed data from acquiring unjustified measurement authority. Consumers must handle both types explicitly, which adds code but preserves uncertainty and scientific meaning.
