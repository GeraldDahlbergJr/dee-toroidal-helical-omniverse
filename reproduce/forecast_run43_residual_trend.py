#!/usr/bin/env python3
"""Diagnostic-only forecast of numerical residual convergence.

Reads a CSV with iteration and positive residual columns. No solver imports,
physical model changes, or automatic changes to acceptance criteria.
Forecasts are heuristic, not guarantees of convergence or failure.
"""
import argparse
import csv
import json
import math
import statistics


def forecast(values, target=1e-10, window=5, horizon=20):
    """Estimate log-residual trend using a sliding linear least-squares fit."""
    if target <= 0 or not math.isfinite(target):
        raise ValueError("target must be finite and positive")
    if window < 3 or horizon < 1:
        raise ValueError("window >= 3 and horizon >= 1 required")
    if len(values) < window:
        return {"status": "insufficient_history", "samples": len(values)}
    if any(v <= 0 or not math.isfinite(v) for v in values):
        raise ValueError("all residuals must be finite and positive")
    recent = values[-window:]
    ys = [math.log(v) for v in recent]
    xs = list(range(window))
    xm = statistics.mean(xs)
    ym = statistics.mean(ys)
    denom = sum((x-xm)**2 for x in xs)
    slope = sum((x-xm)*(y-ym) for x,y in zip(xs,ys))/denom
    fitted = [ym+slope*(x-xm) for x in xs]
    scatter = math.sqrt(sum((y-f)**2 for y,f in zip(ys,fitted))/max(1,window-2))
    ratio = math.exp(slope) if slope < 700 else float("inf")
    latest = recent[-1]
    if latest <= target:
        status, remaining = "target_met", 0
    elif slope >= 0:
        status, remaining = "stagnating_or_growing", None
    else:
        remaining = math.log(target/latest)/slope
        status = "target_unlikely_within_horizon" if remaining > horizon else "target_plausible_within_horizon"
    # Flag noisy fits: a trend prediction is unreliable if variation
    # overwhelms the expected change over the observed window.
    reliable = scatter <= max(abs(slope)*(window-1), 1e-12)
    if not reliable and status != "target_met":
        status = "uncertain_noisy_trend"
    return {"status":status, "samples":len(values), "window":window,
            "latest_residual":latest, "target":target,
            "geometric_reduction_factor":ratio,
            "log_slope_per_iteration":slope,
            "log_fit_scatter":scatter,
            "estimated_iterations_remaining":remaining,
            "forecast_horizon":horizon, "trend_reliable":reliable,
            "note":"Forecast is observational only; always check actual unpreconditioned residual and original gates."}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv_file", help="CSV containing sequential residual history")
    p.add_argument("--column", default="true_relative_residual")
    p.add_argument("--target", type=float, default=1e-10)
    p.add_argument("--window", type=int, default=5)
    p.add_argument("--horizon", type=int, default=20)
    args=p.parse_args()
    with open(args.csv_file, newline="", encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if args.column not in (reader.fieldnames or []):
            p.error("missing column: "+args.column)
        values=[float(row[args.column]) for row in reader if row[args.column].strip()]
    print(json.dumps(forecast(values,args.target,args.window,args.horizon),
                     indent=2,allow_nan=False))


if __name__=="__main__":
    main()
