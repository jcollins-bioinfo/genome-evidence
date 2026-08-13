from genome_evidence.phasing_imputation.masked import masked_metrics, select_masked_markers


def test_masking_is_deterministic_and_reports_denominators() -> None:
    rows = [{"variant_id": f"synthetic-{i}"} for i in range(20)]
    assert select_masked_markers(rows, 0.2, 7) == select_masked_markers(rows, 0.2, 7)
    truth = [{"variant_id": "v1", "reference": "A", "alleles": ["A", "G"]}]
    qc = masked_metrics(truth, [{"variant_id": "v1", "gt": ["G", "A"], "ds": 1.0}])
    assert qc["hard_gt_concordance"] == {"numerator": 1, "denominator": 1}
    assert qc["returned"] == {"numerator": 1, "denominator": 1}
    assert qc["genotype_probability_metrics"] == "not_evaluable"
