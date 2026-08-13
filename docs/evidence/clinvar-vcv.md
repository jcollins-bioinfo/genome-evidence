# ClinVar VCV XML ingestion

Obtain a fixed monthly `ClinVarVariationRelease` VCV XML snapshot using the procedures in the [official ClinVar downloads documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/), then supply that unchanged local `.xml` or `.xml.gz` file. Downloading and network access are deliberately outside this pipeline.

```bash
uv run genome-evidence evidence ingest-clinvar --input /references/fixed.xml.gz --output /private/evidence
uv run genome-evidence evidence link --normalization-run /private/m2 --evidence-run /private/evidence --output /private/annotation
```

The adapter streams `ReleaseSet`/`ClinVarVariationRelease` roots and `VariationArchive` records. It supports simple VCF-style `SequenceLocation` alleles; germline, somatic clinical-impact, and oncogenicity classifications; and submission-level `ClinicalAssertion` data. It preserves exact classification/review/status terms, accessions and versions, submitter data, dates, conditions, citations, method text, and counts of `ObservedIn`/evidence structures. The latter subtrees are **counted, not fully modeled**.

Haplotypes, genotypes, compound records, incomplete locations, non-GRCh38 alleles, and records requiring trimming, liftover, reverse complementation, or repair are not linked. Unknown terms are retained. Missing release identity may be supplied with `--release-identity`; contradictory overrides are rejected. Release date must be present in XML.

Evidence runs contain source metadata, variant/assertion/relationship/condition Parquet tables, aggregate QC/report files, and a checksummed manifest. Annotation runs contain metadata, `variant_evidence_links.parquet`, aggregate QC/report files, and a checksummed manifest. Linking validates upstream manifests and every declared artifact. Outputs are atomic and never overwritten.

Run notebooks with `uv run pytest --nbmake notebooks`; CI executes the same independent, offline, synthetic notebooks.
