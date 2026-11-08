"""Fail closed against frozen 10% and order-1.5 criteria, using rational checks."""
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
data = json.loads((HERE / "error_budget.json").read_text())
scores = {(r["n"], r["region"], r["component"]): r for r in data["scores"]}
assert len(scores) == 60
established = {}
failed = []
for key, row in scores.items():
    allowance = F(row["archived_score_binary64"]) / 10
    assert F(row["max_total_error_exact"]) == allowance
    terms = row["error_budget_absolute_upper"]
    total = row["total_error_upper"]
    ready = all(x is not None for x in terms.values()) and total is not None
    ready &= row["verified_solution_neighborhood"] is not None
    ready &= row["certified_residual_rms_enclosure"] is not None
    if ready:
        assert all(F(x) >= 0 for x in terms.values())
        assert sum(map(F, terms.values()), F(0)) <= F(total)
        if F(total) > allowance:
            failed.append(list(key))
    established[key] = bool(ready and F(total) <= allowance)
conditional = []
certified_pairs = 0
for (n, region, comp), coarse in scores.items():
    if n == 32:
        continue
    fine = scores[2*n, region, comp]
    ratio = 9 * F(coarse["archived_score_binary64"]) / (11 * F(fine["archived_score_binary64"]))
    conditional.append(ratio > 0 and ratio**2 >= 8)
    if established[n, region, comp] and established[2*n, region, comp]:
        actual_ratio = (F(coarse["archived_score_binary64"])-F(coarse["total_error_upper"]))/(F(fine["archived_score_binary64"])+F(fine["total_error_upper"]))
        if actual_ratio > 0 and actual_ratio**2 >= 8:
            certified_pairs += 1
        else:
            failed.append([n, 2*n, region, comp])
assert len(conditional) == 40
# An upper bound exceeding an allowance cannot prove the true error fails.
verdict = "PASS" if all(established.values()) and certified_pairs == 40 else "NOT YET ESTABLISHED"
print(json.dumps({"verdict": verdict, "certified_total_error_scores": sum(established.values()), "score_count": 60, "certified_convergence_comparisons": certified_pairs, "comparison_count": 40, "conditional_10_percent_comparisons_pass": sum(conditional), "proposed_bounds_not_sufficient_for_acceptance": failed}, indent=2))
