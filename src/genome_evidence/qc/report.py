from genome_evidence.qc.models import AssayQCSummary

BOUNDARY = (
    "This report characterizes the supplied genotype file and its source observations. "
    "It does not validate alleles against a reference genome and does not provide medical "
    "or biological interpretation."
)


def render_qc_report(summary: AssayQCSummary) -> str:
    rate = "not applicable" if summary.call_rate is None else f"{summary.call_rate:.6f}"
    return f"""# 23andMe source-file assay QC

> **{BOUNDARY}**

## File and calls

- Source SHA-256: `{summary.source_sha256}`
- Build: `{summary.declared_or_resolved_build}` ({summary.build_provenance.value})
- Parsed / malformed records: {summary.parsed_record_count} / {summary.malformed_record_count}
- Called / no-call records: {summary.called_record_count} / {summary.no_call_record_count}
- Call rate: {rate} (called / all structurally valid parsed records)

## Structural findings

- Duplicate marker occurrences: {summary.duplicate_marker_id_count}
- Duplicate coordinate occurrences: {summary.duplicate_coordinate_count}
- Exact duplicate occurrences: {summary.exact_duplicate_record_count}
- Conflicting duplicate marker occurrences: {summary.conflicting_duplicate_marker_count}
- Unrecognized chromosome-token records: {summary.unrecognized_chromosome_token_count}
- Out-of-order records: {summary.out_of_order_record_count}
"""
