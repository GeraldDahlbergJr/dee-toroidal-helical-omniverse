"""Audit prescribed static radial DEE-action profiles before gravity solves.

These illustrative Gaussian profiles generally DO NOT solve the matter field
equations. Nonzero conservation and field residuals block their use as a
self-consistent stationary harmonic-gauge Einstein source.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from dee_matter import null_contraction, stress_energy


def run(config: dict, output_dir: Path) -> dict:
    R, width = float(config["outer_radius"]), float(config["width"])
    points = int(config["grid_points"])
    if not (np.isfinite(R) and np.isfinite(width) and R > 0 and width > 0 and points >= 5):
        raise ValueError("outer_radius and width must be positive; grid_points >= 5")
    params = {key: float(config[key]) for key in ("lambda_v", "v", "a", "b", "c")}
    amplitudes = np.array([float(config[key]) for key in ("rho_amplitude", "theta_amplitude", "psi_amplitude")])
    if not np.all(np.isfinite(amplitudes)):
        raise ValueError("profile amplitudes must be finite")

    radius = np.linspace(0.0, R, points)
    gaussian = np.exp(-radius**2 / (2 * width**2))
    derivative = -radius / width**2 * gaussian
    rho = params["v"] + amplitudes[0] * gaussian
    field_derivatives = np.outer(derivative, amplitudes)
    gradients = [np.zeros((points, 4)) for _ in range(3)]
    for index in range(3):
        gradients[index][:, 1] = field_derivatives[:, index]
    tensor = stress_energy(rho, *gradients, **params)
    radial_nec = null_contraction(tensor, np.array([1.0, 0.0, 0.0]))
    azimuthal_nec = null_contraction(tensor, np.array([0.0, 1.0, 0.0]))
    # Flat cylindrical covariant divergence for static diagonal stress:
    # ∇_a T^a_r = ∂_r T_rr + (T_rr-T_phiphi)/r.
    conservation = np.zeros(points)
    conservation[1:] = (
        np.gradient(tensor[:, 1, 1], radius, edge_order=2)[1:]
        + (tensor[1:, 1, 1] - tensor[1:, 2, 2]) / radius[1:]
    )

    k_theta = 1 + params["a"] * rho**2
    k_psi = 1 + params["b"] * rho**2
    mixing = params["c"] * rho**2
    theta_prime, psi_prime = field_derivatives[:, 1], field_derivatives[:, 2]
    currents = [k_theta * theta_prime + mixing * psi_prime,
                k_psi * psi_prime + mixing * theta_prime]
    transport_residuals = []
    for current in currents:
        residual = np.zeros(points)
        residual[1:] = np.gradient(radius * current, radius, edge_order=2)[1:] / radius[1:]
        transport_residuals.append(residual)
    rho_prime = field_derivatives[:, 0]
    radial_residual = np.zeros(points)
    radial_residual[1:] = (
        np.gradient(radius * rho_prime, radius, edge_order=2)[1:] / radius[1:]
        - params["lambda_v"] * rho[1:] * (rho[1:]**2 - params["v"]**2)
        - rho[1:] * (params["a"] * theta_prime[1:]**2
                      + params["b"] * psi_prime[1:]**2
                      + 2 * params["c"] * theta_prime[1:] * psi_prime[1:])
    )
    interior = slice(2, -2)
    maxima = {
        "conservation": float(np.max(np.abs(conservation[interior]))),
        "rho_equation": float(np.max(np.abs(radial_residual[interior]))),
        "theta_equation": float(np.max(np.abs(transport_residuals[0][interior]))),
        "psi_equation": float(np.max(np.abs(transport_residuals[1][interior]))),
    }
    summary = {
        "status": "diagnostic prescribed profile; not a self-consistent matter solution",
        "config": config,
        "minimum_radial_nec": float(np.min(radial_nec)),
        "minimum_azimuthal_nec": float(np.min(azimuthal_nec)),
        "max_abs_interior_residuals": maxima,
        "warning": "Nonzero field/conservation residuals preclude using this profile as a stationary harmonic-gauge Einstein source.",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "matter_profile.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["radius", "rho", "T_tt", "T_rr", "T_phiphi", "radial_NEC", "azimuthal_NEC", "conservation_residual"])
        writer.writerows(zip(radius, rho, tensor[:, 0, 0], tensor[:, 1, 1], tensor[:, 2, 2], radial_nec, azimuthal_nec, conservation))
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).parent / "examples" / "matter_profile_config.json")
    parser.add_argument("--output", type=Path, default=Path("matter_profile_output"))
    args = parser.parse_args()
    print(json.dumps(run(json.loads(args.config.read_text()), args.output), indent=2))
