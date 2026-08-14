"""Exact-compatible aggregate reference-score contextualization."""

import bisect
import math
from typing import Any


def compare_reference(
    score: float, reference: dict[str, Any], compatibility: dict[str, str]
) -> dict[str, Any]:
    required = ("pgs_id", "model_version", "assembly", "matching_pipeline", "missingness_policy")
    if any(reference.get(k) != compatibility.get(k) for k in required):
        return {"status": "not_evaluable", "reason": "REFERENCE_IDENTITY_INCOMPATIBLE"}
    values = reference.get("scores")
    if (
        not isinstance(values, list)
        or len(values) < 2
        or not all(isinstance(x, (int, float)) and math.isfinite(x) for x in values)
    ):
        return {"status": "not_evaluable", "reason": "REFERENCE_DISTRIBUTION_INVALID"}
    ordered = sorted(float(x) for x in values)
    mean = math.fsum(ordered) / len(ordered)
    variance = math.fsum((x - mean) ** 2 for x in ordered) / len(ordered)
    if variance <= 0:
        return {"status": "not_evaluable", "reason": "REFERENCE_ZERO_VARIANCE"}
    left, right = bisect.bisect_left(ordered, score), bisect.bisect_right(ordered, score)
    return {
        "status": "evaluable",
        "z_score": (score - mean) / math.sqrt(variance),
        "percentile": 100 * ((left + right) / 2) / len(ordered),
        "sample_count": len(ordered),
        "tail": "lower_or_equal_midrank",
        "tie_method": "midrank",
    }
