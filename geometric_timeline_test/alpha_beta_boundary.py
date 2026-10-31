#!/usr/bin/env python3
"""Isolated alpha/beta phase coupling: no Einstein solver or Run 43 field inputs.

Three spatially oriented, 120-degree-offset phase components form alpha.
Beta is a single linear oscillation. A radial C2 quintic switch compares
sharp versus smooth coupling. Measured orders are derivative errors only.
"""
from pathlib import Path
import json
import math
import numpy as np

OUT = Path(__file__).resolve().parent
R0, WIDTH, TIME = .15, .06, .7
NS = (101, 201, 401, 801)


def alpha(t):
    # Balanced rotating 2D vector from three 120-degree-separated components.
    offsets = 2 * np.pi * np.arange(3) / 3
    directions = np.column_stack((np.cos(offsets), np.sin(offsets)))
    return np.sum(np.cos(t - offsets)[:, None] * directions, axis=0)


def beta(t):
    # Single-phase excitation along a fixed spatial axis.
    return np.array([np.cos(t), 0.])


def switch(r, smooth):
    if not smooth:
        return (r >= R0).astype(float), None, None
    u = np.clip((r - (R0 - WIDTH / 2)) / WIDTH, 0., 1.)
    interior = (u > 0) & (u < 1)
    if smooth == "C3":
        s = 35*u**4 - 84*u**5 + 70*u**6 - 20*u**7
        ds = np.where(interior, (140*u**3 - 420*u**4 + 420*u**5 - 140*u**6) / WIDTH, 0.)
        d2s = np.where(interior, (420*u**2 - 1680*u**3 + 2100*u**4 - 840*u**5) / WIDTH**2, 0.)
    else:
        s = 6 * u**5 - 15 * u**4 + 10 * u**3
        ds = np.where(interior, (30 * u**4 - 60 * u**3 + 30 * u**2) / WIDTH, 0.)
        d2s = np.where(interior, (120 * u**3 - 180 * u**2 + 60 * u) / WIDTH**2, 0.)
    return s, ds, d2s


def main():
    av, bv = alpha(TIME), beta(TIME)
    assert np.linalg.norm(av) > 0. and np.all(np.isfinite(av))
    assert np.linalg.norm(bv) > 0. and np.all(np.isfinite(bv))
    rows = []
    for smooth in (False, True, "C3"):
        for n in NS:
            r = np.linspace(.05, .25, n)
            h = r[1] - r[0]
            s, ds, d2s = switch(r, smooth)
            v = av[None, :] + s[:, None] * (bv - av)[None, :]
            first = (v[2:] - v[:-2]) / (2*h)
            second = (v[2:] - 2*v[1:-1] + v[:-2]) / h**2
            row = dict(profile=("smooth_C3" if smooth == "C3" else "smooth_C2") if smooth else "abrupt", points=n,
                       h=float(h), peak_first=float(np.max(np.linalg.norm(first, axis=1))),
                       peak_second=float(np.max(np.linalg.norm(second, axis=1))))
            if smooth:
                exact_first = ds[1:-1, None] * (bv-av)[None, :]
                exact_second = d2s[1:-1, None] * (bv-av)[None, :]
                row["first_error_rms"] = float(np.sqrt(np.mean(np.sum((first-exact_first)**2, axis=1))))
                row["second_error_rms"] = float(np.sqrt(np.mean(np.sum((second-exact_second)**2, axis=1))))
            rows.append(row)
    smooth_rows = [x for x in rows if x["profile"] == "smooth_C2"]
    c3_rows = [x for x in rows if x["profile"] == "smooth_C3"]
    for i in range(1, len(c3_rows)):
        for key in ("first", "second"):
            err = key + "_error_rms"
            c3_rows[i][key + "_observed_p"] = math.log2(c3_rows[i-1][err]/c3_rows[i][err])
    abrupt_rows = [x for x in rows if x["profile"] == "abrupt"]
    for i in range(1, len(smooth_rows)):
        for key in ("first", "second"):
            err = key + "_error_rms"
            smooth_rows[i][key + "_observed_p"] = math.log2(smooth_rows[i-1][err]/smooth_rows[i][err])
    abrupt_peak_ratio = abrupt_rows[-1]["peak_first"] / abrupt_rows[-2]["peak_first"]
    verdict = dict(
        alpha_is_rotating=bool(np.linalg.norm(alpha(0)-alpha(np.pi/2)) > .1),
        beta_is_single_axis=bool(beta(TIME)[1] == 0),
        abrupt_derivative_spike_grows=bool(abrupt_peak_ratio > 1.8),
        smooth_first_derivative_order_exceeds_1_5=bool(all(x.get("first_observed_p", 2) > 1.5 for x in smooth_rows[1:])),
        C2_second_derivative_below_1_5=bool(all(x["second_observed_p"] < 1.5 for x in smooth_rows[1:])),
        C3_second_derivative_above_1_5=bool(all(x["second_observed_p"] > 1.5 for x in c3_rows[1:])))
    record = dict(scope="Synthetic phase-coupling derivative test only; NOT a Run 43 convergence pass",
                  physical_interpretation="Prescribed excitation; no stress energy, Einstein constraints, or claimed physical cable model",
                  parameters=dict(r0=R0, width=WIDTH, t=TIME, grid_points=NS),
                  rows=rows, checks=verdict, all_checks_passed=all(verdict.values()))
    (OUT / "alpha_beta_boundary_results.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))
    return 0 if all(verdict.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
