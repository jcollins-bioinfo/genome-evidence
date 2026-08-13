"""Deterministic internal masked-marker consistency assessment."""

import hashlib
import math
from typing import Any


def select_masked_markers(
    rows: list[dict[str, Any]], fraction: float, seed: int
) -> tuple[str, ...]:
    if not 0 <= fraction < 1:
        raise ValueError("masked fraction must be in [0, 1)")
    ranked = sorted(
        (hashlib.sha256(f"{seed}:{r['variant_id']}".encode()).hexdigest(), str(r["variant_id"]))
        for r in rows
    )
    count = int(len(ranked) * fraction)
    return tuple(x[1] for x in ranked[:count])


def masked_metrics(truth: list[dict[str, Any]], inferred: list[dict[str, Any]]) -> dict[str, Any]:
    """Exact denominators; unavailable metrics are explicitly not_evaluable."""
    by_id = {str(x["variant_id"]): x for x in inferred}
    returned = [x for x in truth if str(x["variant_id"]) in by_id]
    hard = [x for x in returned if by_id[str(x["variant_id"])].get("gt") is not None]
    correct = sum(sorted(x["alleles"]) == sorted(by_id[str(x["variant_id"])]["gt"]) for x in hard)
    dosage = [
        (sum(a != x["reference"] for a in x["alleles"]), by_id[str(x["variant_id"])].get("ds"))
        for x in returned
    ]
    dosage = [(a, float(b)) for a, b in dosage if b is not None]
    errors = [b - a for a, b in dosage]
    return {
        "assessment": "internal_masked_marker_consistency",
        "returned": {"numerator": len(returned), "denominator": len(truth)},
        "hard_gt_concordance": (
            {"numerator": correct, "denominator": len(hard)} if hard else "not_evaluable"
        ),
        "dosage_mae": sum(abs(x) for x in errors) / len(errors) if errors else "not_evaluable",
        "dosage_rmse": math.sqrt(sum(x * x for x in errors) / len(errors))
        if errors
        else "not_evaluable",
        "genotype_probability_metrics": "not_evaluable",
    }
