"""Reproducible benchmark for the stationary cylindrical radial solver.

The compact polynomial source here is chosen to have a closed-form solution.
It is not derived from the DEE matter action and does not reproduce archived
manuscript CSVs or test a nonlinear Einstein--matter solution.
"""

from __future__ import annotations

import argparse
import csv
import json
import warnings
from pathlib import Path

import numpy as np

from linearized_einstein import (
    LinearizedCylindricalApproximationWarning,
    solve_cylindrical_linearized_field,
)


def compact_source(radius: np.ndarray, outer_radius: float, amplitude: float) -> np.ndarray:
    """Smooth, nonnegative source A(1-(r/R)^2)^2 on the sampled disk."""
    x = radius / outer_radius
    return amplitude * (1.0 - x * x) ** 2


def exact_field(radius: np.ndarray, outer_radius: float, amplitude: float, g: float) -> np.ndarray:
    """Exact solution of (r h')'/r=-16πG T with h(R)=0 and regular axis."""
    r, R = radius, outer_radius
    primitive = (
        (r**2 - R**2) / 4.0
        - (r**4 - R**4) / (8.0 * R**2)
        + (r**6 - R**6) / (36.0 * R**4)
    )
    return -16.0 * np.pi * g * amplitude * primitive


def solve(points: int, outer_radius: float, amplitude: float, g: float) -> dict:
    radius = np.linspace(0.0, outer_radius, points)
    source = compact_source(radius, outer_radius, amplitude)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", LinearizedCylindricalApproximationWarning)
        field = solve_cylindrical_linearized_field(radius, source, gravitational_constant=g)
    exact = exact_field(radius, outer_radius, amplitude, g)
    return {
        "radius": radius,
        "source": source,
        "field": field,
        "exact": exact,
        "max_absolute_error": float(np.max(np.abs(field - exact))),
    }


def run(config: dict, output_dir: Path) -> dict:
    R = float(config["outer_radius"])
    A = float(config["source_amplitude"])
    g = float(config["gravitational_constant"])
    epsilon = float(config["source_scale"])
    grids = [int(n) for n in config["grid_points"]]
    if not (np.isfinite(R) and R > 0 and np.isfinite(A) and A > 0
            and np.isfinite(g) and g > 0 and np.isfinite(epsilon) and epsilon > 0):
        raise ValueError("radius, amplitude, G, and source scale must be finite and positive")
    if len(grids) < 2 or any(n < 3 for n in grids) or any(b <= a for a, b in zip(grids, grids[1:])):
        raise ValueError("grid_points must contain at least two strictly increasing sizes >= 3")

    runs = [solve(n, R, A, g) for n in grids]
    errors = [item["max_absolute_error"] for item in runs]
    finest = runs[-1]
    scaled_source = epsilon * finest["source"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", LinearizedCylindricalApproximationWarning)
        scaled_field = solve_cylindrical_linearized_field(
            finest["radius"], scaled_source, gravitational_constant=g
        )
    scaling_error = float(np.max(np.abs(scaled_field - epsilon * finest["field"])))
    summary = {
        "status": "solver benchmark only; source is illustrative, not DEE-derived",
        "equation": "(r h')'/r = -16*pi*G*T(r); h(R)=0; regular axis",
        "source": "T(r) = A*(1-(r/R)^2)^2, 0<=r<=R",
        "config": config,
        "max_absolute_error_by_grid": dict(zip(map(str, grids), errors)),
        "error_ratio_coarse_to_fine": [a / b for a, b in zip(errors, errors[1:])],
        "largest_scaled_field": float(np.max(np.abs(scaled_field))),
        "source_scaling_absolute_difference": scaling_error,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "radial_profile.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["radius", "source", "computed_h", "exact_h", "scaled_h"])
        writer.writerows(zip(finest["radius"], finest["source"], finest["field"], finest["exact"], scaled_field))
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).parent / "examples" / "weak_field_config.json")
    parser.add_argument("--output", type=Path, default=Path("weak_field_output"))
    args = parser.parse_args()
    print(json.dumps(run(json.loads(args.config.read_text()), args.output), indent=2))
