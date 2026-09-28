"""Stationary axisymmetric linearized-gravity twist-sector solver.

For a cylindrical covariant perturbation h_{phi z}=H(r), the flat-background
harmonic-gauge tensor equation reduces, under stationarity, axisymmetry and
z-independence, to

    H'' - H'/r = -16 pi G T_{phi z}.

This module implements that restricted equation only. It is not a nonlinear GR
solver and does not establish physical validity of a complete DEE model.
"""

from __future__ import annotations

import numpy as np


def solve_tensor_twist(
    radius: np.ndarray,
    stress_phi_z: np.ndarray,
    gravitational_constant: float = 1.0,
    outer_boundary: float = 0.0,
) -> np.ndarray:
    """Solve the tensor twist BVP with H(0)=0 and H(R)=outer_boundary.

    The grid must be uniform and start at r=0. H(0)=0 is the coordinate-basis
    regularity condition implied by smooth Cartesian components for the tested
    axisymmetric sector. Interior derivatives use centered second-order finite
    differences.
    """
    r = np.asarray(radius, dtype=float)
    stress = np.asarray(stress_phi_z, dtype=float)
    if r.ndim != 1 or r.size < 3:
        raise ValueError("radius must be one-dimensional with at least 3 points")
    if stress.shape != r.shape:
        raise ValueError("stress_phi_z must have the same shape as radius")
    if not np.all(np.isfinite(r)) or not np.all(np.isfinite(stress)):
        raise ValueError("radius and stress_phi_z must be finite")
    if r[0] != 0.0 or np.any(np.diff(r) <= 0.0):
        raise ValueError("radius must start at 0 and be strictly increasing")
    drs = np.diff(r)
    if not np.allclose(drs, drs[0], rtol=1e-12, atol=1e-15):
        raise ValueError("radius grid must be uniform")
    if not np.isfinite(gravitational_constant) or gravitational_constant < 0.0:
        raise ValueError("gravitational_constant must be finite and non-negative")
    if not np.isfinite(outer_boundary):
        raise ValueError("outer_boundary must be finite")

    n = r.size
    dr = drs[0]
    matrix = np.zeros((n, n), dtype=float)
    rhs = -16.0 * np.pi * gravitational_constant * stress.copy()

    matrix[0, 0] = 1.0
    rhs[0] = 0.0
    matrix[-1, -1] = 1.0
    rhs[-1] = outer_boundary

    for i in range(1, n - 1):
        ri = r[i]
        matrix[i, i - 1] = 1.0 / dr**2 + 1.0 / (2.0 * dr * ri)
        matrix[i, i] = -2.0 / dr**2
        matrix[i, i + 1] = 1.0 / dr**2 - 1.0 / (2.0 * dr * ri)

    return np.linalg.solve(matrix, rhs)
