# M1 assay-level QC

M1 QC is descriptive source-file QC, not biological or clinical quality classification. It has no global PASS/FAIL and no low-call-rate threshold.

`eligible_parsed_record_count` is every structurally valid four-column record with a positive integer source position. Each eligible record is exactly `CALLED` or `NO_CALL`; malformed rows are excluded and separately counted. Therefore:

```text
call_rate = called_record_count / eligible_parsed_record_count
```

The rate is null when the denominator is zero. Per-source-chromosome rates use the analogous per-token denominator. Chromosome recognition is lexical only (`1`–`22`, `X`, `Y`, `XY`, `M`, `MT`); other tokens remain unchanged and receive findings. Lexical homozygous/heterozygous labels apply only to two-character A/C/G/T tokens and do not assert biological ploidy or zygosity.

Duplicate metrics count occurrences after the first within each duplicate group. Exact repeats, duplicate marker IDs, duplicate coordinates, and marker conflicts are reported independently and retained. Ordering diagnostics count a position decrease relative to the prior record with the same exact source chromosome token. Findings carry codes, severity, and bounded line references rather than genotype dumps.
