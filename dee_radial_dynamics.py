"""Time-dependent radial evolution of the manuscript's three matter fields.

Finite-volume radial fluxes and a semidiscrete Hamiltonian use the polynomial
matter action on a fixed flat cylindrical background. This is a matter-only
model: it neither solves Einstein's equations nor implements angular winding.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from dee_matter import null_contraction, stress_energy


def _coefficients(rho: np.ndarray, config: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return 1 + config["a"] * rho**2, 1 + config["b"] * rho**2, config["c"] * rho**2


def _kinematics(state: np.ndarray, config: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rho, _, _, _, ptheta, ppsi = state
    kt, kp, mix = _coefficients(rho, config)
    determinant = kt * kp - mix**2
    if np.any(kt <= 0) or np.any(determinant <= 0):
        raise ValueError("transport kinetic matrix lost positive definiteness")
    return (kp * ptheta - mix * ppsi) / determinant, (kt * ppsi - mix * ptheta) / determinant, determinant


def rhs(state: np.ndarray, config: dict) -> np.ndarray:
    """Hamilton equations for cell-centered (rho,theta,psi,p_rho,p_theta,p_psi)."""
    rho, theta, psi, prho, _, _ = state
    n = rho.size
    dr = config["outer_radius"] / n
    radii = (np.arange(n) + 0.5) * dr
    faces = np.arange(n + 1) * dr
    vtheta, vpsi, _ = _kinematics(state, config)

    drho, dtheta, dpsi = ((field[1:] - field[:-1]) / dr for field in (rho, theta, psi))
    rhoface = (rho[1:] + rho[:-1]) / 2
    kt, kp, mix = _coefficients(rhoface, config)
    fluxes = np.zeros((3, n + 1))  # zero flux at the axis and reflecting outer face
    fluxes[0, 1:-1] = drho
    fluxes[1, 1:-1] = kt * dtheta + mix * dpsi
    fluxes[2, 1:-1] = kp * dpsi + mix * dtheta
    divergence = (faces[1:] * fluxes[:, 1:] - faces[:-1] * fluxes[:, :-1]) / (radii * dr)

    # Variation of the face-centered K(rho) gradient energy. A face shares
    # its dependence on rho equally between its neighboring cells.
    d_k_quadratic = 2 * rhoface * (
        config["a"] * dtheta**2 + config["b"] * dpsi**2
        + 2 * config["c"] * dtheta * dpsi
    )
    face_term = faces[1:-1] * d_k_quadratic
    radial_coupling = np.zeros(n)
    radial_coupling[:-1] += face_term
    radial_coupling[1:] += face_term
    radial_coupling /= 4 * radii
    potential_derivative = config["lambda_v"] * rho * (rho**2 - config["v"]**2)
    temporal_coupling = rho * (
        config["a"] * vtheta**2 + config["b"] * vpsi**2
        + 2 * config["c"] * vtheta * vpsi
    )
    return np.stack((
        prho, vtheta, vpsi,
        divergence[0] - potential_derivative + temporal_coupling - radial_coupling,
        divergence[1], divergence[2],
    ))


def total_energy(state: np.ndarray, config: dict) -> float:
    """Discrete matter energy, including the same face coupling used by rhs."""
    rho, theta, psi, prho, ptheta, ppsi = state
    n = rho.size
    dr = config["outer_radius"] / n
    radii = (np.arange(n) + 0.5) * dr
    faces = np.arange(1, n) * dr
    vtheta, vpsi, _ = _kinematics(state, config)
    potential = config["lambda_v"] / 4 * (rho**2 - config["v"]**2) ** 2
    temporal = 0.5 * (prho**2 + ptheta * vtheta + ppsi * vpsi) + potential
    drho, dtheta, dpsi = ((field[1:] - field[:-1]) / dr for field in (rho, theta, psi))
    kt, kp, mix = _coefficients((rho[1:] + rho[:-1]) / 2, config)
    spatial = 0.5 * (drho**2 + kt * dtheta**2 + kp * dpsi**2 + 2 * mix * dtheta * dpsi)
    return float(2 * np.pi * dr * (np.dot(radii, temporal) + np.dot(faces, spatial)))


def initial_state(config: dict) -> np.ndarray:
    n = config["grid_points"]
    r = (np.arange(n) + 0.5) * config["outer_radius"] / n
    envelope = np.exp(-r**2 / (2 * config["width"]**2))
    state = np.zeros((6, n))
    state[0] = config["v"] + config["rho_amplitude"] * envelope
    state[1] = config["theta_amplitude"] * envelope
    state[2] = config["psi_amplitude"] * envelope
    return state


def _snapshot(state: np.ndarray, time: float, config: dict) -> tuple[list, dict]:
    n = state.shape[1]
    r = (np.arange(n) + 0.5) * config["outer_radius"] / n
    dr = config["outer_radius"] / n
    vtheta, vpsi, determinant = _kinematics(state, config)
    gradients = [np.zeros((n, 4)) for _ in range(3)]
    gradients[0][:, 0] = state[3]
    gradients[1][:, 0] = vtheta
    gradients[2][:, 0] = vpsi
    for index in range(3):
        gradients[index][:, 1] = np.gradient(state[index], dr, edge_order=2)
    params = {key: config[key] for key in ("lambda_v", "v", "a", "b", "c")}
    tensor = stress_energy(state[0], *gradients, **params)
    radial_nec = null_contraction(tensor, np.array([1.0, 0.0, 0.0]))
    azimuthal_nec = null_contraction(tensor, np.array([0.0, 1.0, 0.0]))
    rows = list(zip(np.full(n, time), r, state[0], state[1], state[2],
                    tensor[:, 0, 0], tensor[:, 0, 1], radial_nec, azimuthal_nec))
    metadata = {
        "time": time,
        "energy": total_energy(state, config),
        "minimum_rho": float(np.min(state[0])),
        "minimum_kinetic_determinant": float(np.min(determinant)),
        "minimum_radial_nec": float(np.min(radial_nec)),
        "minimum_azimuthal_nec": float(np.min(azimuthal_nec)),
        "max_energy_density": float(np.max(tensor[:, 0, 0])),
    }
    return rows, metadata


def run(config: dict, output_dir: Path) -> dict:
    config = dict(config)
    for key in ("outer_radius", "width", "lambda_v", "v", "a", "b", "c",
                "rho_amplitude", "theta_amplitude", "psi_amplitude", "time_step"):
        config[key] = float(config[key])
    for key in ("grid_points", "steps", "snapshot_interval"):
        config[key] = int(config[key])
    if not all(np.isfinite(value) for value in config.values()):
        raise ValueError("all configuration values must be finite")
    if (min(config[key] for key in ("outer_radius", "width", "lambda_v", "v", "time_step")) <= 0
            or config["grid_points"] < 5 or config["steps"] < 1 or config["snapshot_interval"] < 1):
        raise ValueError("positive scales and at least five cells are required")
    if config["time_step"] >= config["outer_radius"] / config["grid_points"] / 2:
        raise ValueError("time_step must be less than half the cell width")

    state = initial_state(config)
    rows, snapshots = [], []
    for step in range(config["steps"] + 1):
        if step % config["snapshot_interval"] == 0 or step == config["steps"]:
            profile, info = _snapshot(state, step * config["time_step"], config)
            rows.extend(profile)
            snapshots.append(info)
        if step == config["steps"]:
            break
        h = config["time_step"]
        k1 = rhs(state, config)
        k2 = rhs(state + 0.5 * h * k1, config)
        k3 = rhs(state + 0.5 * h * k2, config)
        k4 = rhs(state + h * k3, config)
        state = state + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    energies = np.array([snapshot["energy"] for snapshot in snapshots])
    summary = {
        "status": "time-dependent matter-only evolution; no gravitational solve",
        "config": config,
        "snapshots": snapshots,
        "max_relative_energy_drift": float(np.max(np.abs(energies - energies[0])) / energies[0]) if energies[0] else 0.0,
        "warning": "The stationary radial Einstein equation cannot use these time-dependent sources as a solved spacetime.",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "dynamic_profiles.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "radius", "rho", "theta", "psi", "T_tt", "T_tr", "radial_NEC", "azimuthal_NEC"])
        writer.writerows(rows)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).parent / "examples" / "dynamic_matter_config.json")
    parser.add_argument("--output", type=Path, default=Path("dynamic_matter_output"))
    args = parser.parse_args()
    print(json.dumps(run(json.loads(args.config.read_text()), args.output), indent=2))
