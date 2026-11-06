"""Exact rational checks on archived scores and certified Schur inputs."""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

source = Path(sys.argv[1])
rows = json.loads((source / "results.json").read_text())["rows"]
blocks = json.loads((source / "conditioning.json").read_text(), parse_float=lambda value: Fraction(float(value)))
assert [row["n"] for row in rows] == [8, 16, 32]
budget = Fraction(1, 10)
inverse_cap = Fraction(133, 100)
out = {"verdict": "NOT YET ESTABLISHED", "scores": [], "pairs": []}
for row, block in zip(rows, blocks, strict=True):
    a = block["scalar_gradient_coercivity_lower"]
    v = block["vector_gradient_coercivity_lower"]
    c = block["scalar_cross_gradient_norm_upper"]
    d = block["vector_cross_gradient_norm_upper"]
    s = block["schur_gradient_norm_upper"]
    scalar = (1 + c / v) / (a * (1 - s))
    vector = (1 + d * scalar) / v
    assert max(scalar, vector) < inverse_cap
    for region, values in row["evaluator_q8"]["residual_rms"].items():
        for i, (component, value) in enumerate(zip(("H", "Mx", "My", "Mz"), values, strict=True)):
            score = Fraction(value)
            assert score > 0
            q6 = Fraction(row["evaluator_q6"]["residual_rms"][region][i])
            out["scores"].append({
                "n": row["n"], "region": region, "component": component,
                "archived_score_binary64": value,
                "max_total_error_exact": str(score * budget),
                "empirical_proxy_relative": float(2 * abs(score - q6) / score),
                "sampled_inverse_upper_approx": float(max(scalar, vector)),
                "certified_sampled_inverse_cap": 1.33,
                "total_error_upper": None, "verdict": "NOT YET ESTABLISHED",
            })
scores = {(x["n"], x["region"], x["component"]): x for x in out["scores"]}
assert len(scores) == 60
for (n, region, component), coarse in scores.items():
    if n == 32:
        continue
    fine = scores[2*n, region, component]
    ratio = (1-budget)*Fraction(coarse["archived_score_binary64"])/(
        (1+budget)*Fraction(fine["archived_score_binary64"]))
    out["pairs"].append({
        "pair": [n, 2*n], "region": region, "component": component,
        "conditional_ratio_exact": str(ratio),
        "conditional_order_display": math.log2(float(ratio)),
        "conditional_order_at_least_1_5": ratio > 0 and ratio**2 >= 8,
        "certified_order": None,
    })
assert len(out["pairs"]) == 40
assert all(p["conditional_order_at_least_1_5"] for p in out["pairs"])
print(json.dumps(out, indent=2))
