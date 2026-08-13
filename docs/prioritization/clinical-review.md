# Clinical review routing

M4 builds a private research review queue from checksum-valid M2 normalization, M3 evidence, and M3 exact-link runs. It runs offline and requires an explicit policy and `germline_constitutional` analysis context.

```bash
genome-evidence prioritize clinical \
  --normalization-run /private/m2-run \
  --evidence-run /private/evidence-run \
  --annotation-run /private/annotation-run \
  --policy references/clinvar-germline-review-policy-v1.json \
  --analysis-context germline_constitutional \
  --output /private/m4-run
```

The shipped `clinvar-germline-review-routing` policy version `1.0.0` routes active germline high-attention source terms to `review_first`; risk, uncertainty, other, unmapped, or source-reported conflict context to `review_next`; benign-only evidence to `not_prioritized`; and somatic-only or inactive evidence to `context_only`. Reference-only and missing canonical genotypes are `not_eligible`. Discordant calls remain unresolved data conflicts.

The output directory is atomic and immutable and contains profiles, candidates, assertion links, structured rationales, exclusions, QC, metadata, a manifest, and a private Markdown report. Every declared artifact is SHA-256 registered. Stable scientific IDs can repeat across identical inputs and policy; runtime run IDs and timestamps intentionally cannot.

Source classifications are not project conclusions. Review bands do not establish disease, probability, clinical validity, actionability, or medical advice. Array data are not comprehensive and absence from the queue is not a negative genetic test.

## Policy matching and deterministic ordering

The default policy matches source terms by Unicode case folding plus collapsed whitespace while
retaining the exact source string in every profile and assertion link. Unknown terms remain
`unmapped_source_term`; unknown record status is active only because policy version `1.0.0` says
so and therefore produces a warning rationale. Replaced and deleted assertions remain visible but
inactive. Within each band, the emitted ordering tuple uses source review level as contextual
metadata, source-reported conflict, missing evaluation date, age measured against the snapshot
release date, and finally the stable profile ID. These fields are ordering inputs, not probabilities
or hidden numeric clinical scores.
