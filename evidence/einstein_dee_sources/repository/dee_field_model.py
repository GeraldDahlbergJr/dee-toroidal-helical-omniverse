"""Minimal finite-equilibrium relations for the DEE model.

The functions in this module implement algebraic relations used for numerical
exploration.  They are not a derivation of a complete physical theory.
"""

from __future__ import annotations

import numpy as np


def _validate_rho_squared(rho_squared: float) -> None:
    """Require a finite, non-negative value for a quantity denoted ``rho^2``."""
    if not np.isfinite(rho_squared) or rho_squared < 0.0:
        raise ValueError("rho_squared must be a finite non-negative value")


def _validate_lambda_v(lambda_v: float) -> None:
    """Require a finite, strictly positive quartic potential coupling."""
    if not np.isfinite(lambda_v) or lambda_v <= 0.0:
        raise ValueError("lambda_v must be a finite strictly positive value")


def kinetic_matrix(a: float, b: float, c: float, rho_squared: float) -> np.ndarray:
    """Return the field-dependent minimal-DEE kinetic matrix at ``rho_squared``.

    The diagonal kinetic coefficients are ``K_Theta = 1 + a*rho_squared``
    and ``K_Psi = 1 + b*rho_squared``; the mixed coefficient is
    ``Lambda = c*rho_squared``.
    """
    _validate_rho_squared(rho_squared)
    return np.array(
        [
            [1.0 + a * rho_squared, c * rho_squared],
            [c * rho_squared, 1.0 + b * rho_squared],
        ],
        dtype=float,
    )


def kinetic_determinant(a: float, b: float, c: float, rho_squared: float) -> float:
    """Return ``det(K)`` at the supplied field value.

    ``det(K) = 1 + (a + b) rho_squared + (a*b - c**2) rho_squared**2``.
    """
    _validate_rho_squared(rho_squared)
    return float(
        1.0
        + (a + b) * rho_squared
        + (a * b - c * c) * rho_squared * rho_squared
    )


def is_positive_definite(a: float, b: float, c: float, rho_squared: float) -> bool:
    """Apply Sylvester's criterion to the field-dependent kinetic matrix."""
    _validate_rho_squared(rho_squared)
    k_theta = 1.0 + a * rho_squared
    return bool(k_theta > 0.0 and kinetic_determinant(a, b, c, rho_squared) > 0.0)

def finite_equilibrium_squared(
    v: float,
    a: float,
    b: float,
    c: float,
    x_theta: float,
    x_psi: float,
    x_theta_psi: float,
    lambda_v: float,
) -> float:
    """Compute ``rho_*^2`` for the minimal finite-equilibrium branch.

    ``rho_*^2 = v^2 - (a X_theta + b X_psi + 2 c X_thetaPsi) / lambda_V``.
    The quartic potential coupling ``lambda_v`` must be finite and strictly
    positive for this minimal branch.
    """
    _validate_lambda_v(lambda_v)
    transport = a * x_theta + b * x_psi + 2.0 * c * x_theta_psi
    return float(v * v - transport / lambda_v)


def finite_equilibrium(*args: float) -> float:
    """Return the non-negative finite-equilibrium amplitude ``rho_*``.

    Raises:
        ValueError: if the requested branch has negative ``rho_*^2``.
    """
    rho_squared = finite_equilibrium_squared(*args)
    if rho_squared < 0.0:
        raise ValueError("finite-equilibrium branch requires rho_*^2 >= 0")
    return float(np.sqrt(rho_squared))


def effective_mass_squared(lambda_v: float, rho_star_squared: float) -> float:
    """Compute the branch curvature ``m_eff^2 = 2 lambda_V rho_*^2``."""
    _validate_lambda_v(lambda_v)
    _validate_rho_squared(rho_star_squared)
    return float(2.0 * lambda_v * rho_star_squared)
