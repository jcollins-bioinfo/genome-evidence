# 23andMe raw-genotype ingestion (M1)

The importable API is:

```python
from pathlib import Path
from genome_evidence.ingest import Ingest23andMeConfig, ParseMode, ingest_23andme

result = ingest_23andme(
    Path("/private/input.txt"),
    Path("/private/run"),
    Ingest23andMeConfig(mode=ParseMode.STRICT, sample_id="pseudonymous-sample"),
)
```

Strict mode is the default for formal processing and stops at the first structurally malformed data line with a privacy-safe line-numbered exception. Lenient mode preserves structured malformed-record diagnostics and continues with valid records. Both modes read bytes once, hash the exact input, tolerate UTF-8 BOM and LF/CRLF transport encoding, preserve source ordering, tokens, comment lines, and line numbers, and perform approximately O(n) parsing/QC. Records are never deduplicated.

A recognized explicit `build` or `assembly` comment is vendor-declared metadata. Without it the build is `UNKNOWN`; coordinates are never used to guess. `genome_build_override` is separately tagged `USER_OVERRIDE` and the declared value remains independently recorded.

## M1 establishes

- what records and comments the file contains;
- source metadata, exact-file SHA-256, and line provenance;
- source-level genotype observations and explicit no-call state;
- structural coherence, duplicate/order diagnostics, and descriptive assay QC.

## M1 does not establish

M1 does not establish biological REF/ALT, canonical variants or genotypes, reference-genome correctness, strand orientation, pathogenicity, disease risk, ancestry, or imputed variants. In particular, `raw source genotype != canonical biological genotype`.

## Outputs and determinism

A new/non-empty-safe private output directory receives `manifest.json`, `source_metadata.json`, `observations.parquet`, `qc_summary.json`, `qc_findings.parquet`, and `qc_report.md`. Inside the Git worktree, outputs are rejected unless Git confirms the path is ignored. The source is never copied. JSON is key-sorted and Parquet retains source ordering. Run IDs and timestamps are execution metadata and therefore vary; aggregate QC and the Markdown report are deterministic for identical source/configuration.
