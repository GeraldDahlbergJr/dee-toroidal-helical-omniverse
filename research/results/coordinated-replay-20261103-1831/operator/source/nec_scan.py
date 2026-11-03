"""Null-energy-condition quadratic-form diagnostics for the minimal DEE model."""

from __future__ import annotations

import numpy as np


def null_energy_quadratic_form(kinetic: np.ndarray, transport: np.ndarray) -> float:
    """Evaluate ``T_kk = u.T @ K @ u`` for one two-component transport vector."""
    matrix = np.asarray(kinetic, dtype=float)
    vector = np.asarray(transport, dtype=float)
    if matrix.shape != (2, 2):
        raise ValueError("kinetic must have shape (2, 2)")
    if vector.shape != (2,):
        raise ValueError("transport must have shape (2,)")
    return float(vector @ matrix @ vector)


def scan_transport_projections(kinetic: np.ndarray, samples: int = 360) -> dict[str, np.ndarray | float | bool]:
    """Sample unit transport directions and report their NEC quadratic forms.

    The scan is diagnostic rather than a replacement for an all-directions
    analytic proof.  It samples ``u=(cos(phi), sin(phi))`` over one full turn.
    """
    if samples < 3:
        raise ValueError("samples must be at least 3")
    matrix = np.asarray(kinetic, dtype=float)
    if matrix.shape != (2, 2):
        raise ValueError("kinetic must have shape (2, 2)")

    angles = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    transports = np.column_stack((np.cos(angles), np.sin(angles)))
    t_kk = np.einsum("ni,ij,nj->n", transports, matrix, transports)
    return {
        "angles": angles,
        "transports": transports,
        "t_kk": t_kk,
        "minimum": float(t_kk.min()),
        "all_nonnegative": bool(np.all(t_kk >= 0.0)),
    }
