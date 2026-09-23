"""Radial solver for a stationary, cylindrical *linearized* gravity model.

This module solves the stated weak-field radial equation only.  It is not a
nonlinear general-relativity solver and cannot establish nonlinear stability,
global non-collapse, or an exact spacetime solution.
"""

from __future__ import annotations

import warnings

import numpy as np


class LinearizedCylindricalApproximationWarning(UserWarning):
    """Signals the restricted physical scope of the radial calculation."""


def solve_cylindrical_linearized_field(
    radius: np.ndarray,
    stress_energy: np.ndarray,
    gravitational_constant: float = 1.0,
    outer_boundary: float | np.ndarray = 0.0,
) -> np.ndarray:
    """Integrate ``(r h')'/r = -16 pi G T(r)`` inward from a finite boundary.

    ``radius`` must begin at zero and be strictly increasing.  Regularity at
    the axis is imposed through ``r h'(r) -> 0`` (and hence ``h'(0)=0``).
    The supplied finite ``outer_boundary`` sets ``h(radius[-1])``.  A trailing
    component shape on ``stress_energy`` is retained, allowing independent
    numerical solves for multiple ``T_mu_nu`` components.
    """
    warnings.warn(
        "This is a linearized stationary cylindrical approximation, not a "
        "nonlinear GR solution.",
        LinearizedCylindricalApproximationWarning,
        stacklevel=2,
    )
    r = np.asarray(radius, dtype=float)
    source_stress = np.asarray(stress_energy, dtype=float)
    boundary = np.asarray(outer_boundary, dtype=float)
    if r.ndim != 1 or r.size < 3:
        raise ValueError("radius must be a one-dimensional grid with at least 3 points")
    if not np.all(np.isfinite(r)):
        raise ValueError("radius must contain only finite values")
    if r[0] != 0.0 or np.any(np.diff(r) <= 0.0):
        raise ValueError("radius must start at 0 and be strictly increasing")
    if source_stress.ndim < 1 or source_stress.shape[0] != r.size:
        raise ValueError("stress_energy first dimension must match radius")
    if not np.all(np.isfinite(source_stress)):
        raise ValueError("stress_energy must contain only finite values")
    if not np.isfinite(gravitational_constant) or gravitational_constant < 0.0:
        raise ValueError("gravitational_constant must be a finite non-negative value")
    if not np.all(np.isfinite(boundary)):
        raise ValueError("outer_boundary must contain only finite values")
    try:
        boundary_components = np.broadcast_to(boundary, source_stress.shape[1:])
    except ValueError as error:
        raise ValueError(
            "outer_boundary must be scalar or broadcastable to stress_energy components"
        ) from error

    source = -16.0 * np.pi * gravitational_constant * source_stress
    integrand = r.reshape((r.size,) + (1,) * (source.ndim - 1)) * source
    radial_flux = np.zeros_like(source)
    increments = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r).reshape(
        (r.size - 1,) + (1,) * (source.ndim - 1)
    )
    radial_flux[1:] = np.cumsum(increments, axis=0)

    derivative = np.zeros_like(source)
    derivative[1:] = radial_flux[1:] / r[1:].reshape((r.size - 1,) + (1,) * (source.ndim - 1))
    field = np.empty_like(source)
    field[-1] = boundary_components
    for index in range(r.size - 2, -1, -1):
        field[index] = field[index + 1] - 0.5 * (derivative[index] + derivative[index + 1]) * (r[index + 1] - r[index])
    return field
