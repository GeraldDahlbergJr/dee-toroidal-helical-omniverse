"""Reproduce a documented weak-field cylindrical solver checkpoint.

This script intentionally tests only the equation implemented by
``solve_cylindrical_linearized_field``. It does not validate nonlinear GR,
DEE as a complete physical theory, or the candidate tensor twist operator.

Run from the repository root:
    python reproduce/weak_field_checkpoint.py
"""

from __future__ import annotations

import json
import pathlib
import sys
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from linearized_einstein import (  # noqa: E402
    LinearizedCylindricalApproximationWarning,
    solve_cylindrical_linearized_field,
)

G = 0.2
T0 = 1.3
ALPHA = 0.7
R_MAX = 2.0
RESOLUTIONS = (41, 81, 161, 321)


def stress_profile(radius: np.ndarray) -> np.ndarray:
    """Manufactured stationary source T(r)=T0(1+alpha r^2)."""
    return T0 * (1.0 + ALPHA * radius**2)


def exact_field(radius: np.ndarray) -> np.ndarray:
    """Exact solution with h(R_MAX)=0 for the manufactured source."""
    scale = -16.0 * np.pi * G * T0
    return scale * (
        (radius**2 - R_MAX**2) / 4.0
        + ALPHA * (radius**4 - R_MAX**4) / 16.0
    )


def observed_order(coarse_error: float, fine_error: float) -> float:
    return float(np.log(coarse_error / fine_error) / np.log(2.0))


def main() -> None:
    rows = []
    previous_error = None
    for points in RESOLUTIONS:
        radius = np.linspace(0.0, R_MAX, points)
        stress = stress_profile(radius)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", LinearizedCylindricalApproximationWarning)
            field = solve_cylindrical_linearized_field(
                radius, stress, gravitational_constant=G, outer_boundary=0.0
            )

        error = float(np.max(np.abs(field - exact_field(radius))))
        row = {
            "points": points,
            "dr": float(radius[1] - radius[0]),
            "max_abs_error": error,
            "outer_boundary_abs_error": float(abs(field[-1])),
            "axis_first_increment": float(field[1] - field[0]),
            "observed_order": None,
        }
        if previous_error is not None:
            row["observed_order"] = observed_order(previous_error, error)
        rows.append(row)
        previous_error = error

    payload = {
        "equation": "(r h')'/r = -16 pi G T(r)",
        "source": "T(r) = T0 (1 + alpha r^2)",
        "parameters": {"G": G, "T0": T0, "alpha": ALPHA, "R_max": R_MAX},
        "boundary_conditions": ["regular axis", "h(R_max)=0"],
        "resolutions": rows,
        "acceptance": {
            "target_observed_order": 2.0,
            "minimum_accepted_order": 1.9,
            "all_refinements_pass": all(
                row["observed_order"] is None or row["observed_order"] >= 1.9
                for row in rows
            ),
        },
    }

    output = ROOT / "reproduce" / "weak_field_checkpoint.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"\nWrote {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
