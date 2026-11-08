"""Build explicit unresolved budgets; certify only denominator bounds.

No PDE solve and no empirical difference is used as an error enclosure.
All proof arithmetic below is exact rational arithmetic on stored IEEE values.
"""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "research/run43_bspline_completion"
PI_LO = F("3.14159265358979323846264338327950288419716939937510")
PI_HI = PI_LO + F(1, 10**50)


def verify_pi_bracket():
    # Alternating arctan series: consecutive partial sums enclose atan(x).
    # Machin's identity gives pi = 16 atan(1/5) - 4 atan(1/239).
    def atan_bounds(x, terms):
        s = sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms)), F(0))
        adjacent = s + (-1)**terms*x**(2*terms+1)/F(2*terms+1)
        return min(s, adjacent), max(s, adjacent)
    a, b = atan_bounds(F(1, 5), 40)
    c, d = atan_bounds(F(1, 239), 12)
    assert PI_LO <= 16*a-4*d <= 16*b-4*c <= PI_HI


def build():
    verify_pi_bracket()
    frozen = json.loads((ROOT / "research/run43_final_gate/audit.json").read_text())
    rows = []
    for n in (8, 16, 32):
        path = SOURCE / f"coefficients_{n}.npz"
        with np.load(path, allow_pickle=False) as z:
            c = z["coefficients"]
            assert int(z["n"]) == n and c.shape == (4, n+2, n, n)
            assert c.dtype == np.float64 and np.isfinite(c).all()
            lo = 1 + min(F(0), F(float(c[0].min())))
            hi = 1 + max(F(0), F(float(c[0].max())))
        assert lo > 0
        for score in frozen["scores"]:
            if score["n"] != n:
                continue
            region = score["region"]
            if region.startswith("first_"):
                r = F(1, 20)+F(1, 5*n) if region == "first_inner" else F(1, 4)-F(1, 5*n)
                radial = r
            else:
                a, b = {"inner": (F(1,20), F(1,10)), "bulk": (F(1,10), F(1,5)), "outer": (F(1,5), F(1,4))}[region]
                radial = (b*b-a*a)/2
            denominator = [4*PI_LO**2*radial*lo**6, 4*PI_HI**2*radial*hi**6]
            rows.append({
                **score,
                "coefficient_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "continuous_weight_denominator_exact_rational_enclosure": list(map(str, denominator)),
                "denominator_upper_over_lower_display": float(denominator[1]/denominator[0]),
                "error_budget_absolute_upper": {
                    "quadrature_truncation": None,
                    "floating_point_derivatives": None,
                    "floating_point_geometry": None,
                    "floating_point_source": None,
                    "floating_point_integration": None,
                    "propagated_exact_minus_gauss6_assembly": None,
                    "propagated_algebraic_solver": None,
                    "nonlinear_remainder": None,
                },
                "certified_residual_rms_enclosure": None,
                "verified_solution_neighborhood": None,
            })
    assert len(rows) == 60
    return {"frozen_checkpoint": "c5580eea5b216052ba5e91c6725e0304a7852a99", "verdict": "NOT YET ESTABLISHED", "denominator_target": "Exact decimal protocol chart and ideal degree-four spline knots; does not enclose binary64 evaluator knot perturbation.", "null_semantics": "Unestablished bound, never zero.", "scores": rows, "convergence_comparisons": frozen["pairs"]}


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))
