"""Flat-background stress tensor of manuscript Eqs. (3), (8), and (9).

Components are covariant in a local orthonormal frame with signature (-,+,+,+).
This algebraic evaluator does not enforce the scalar equations of motion.
"""

from __future__ import annotations

import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def stress_energy(
    rho: np.ndarray,
    grad_rho: np.ndarray,
    grad_theta: np.ndarray,
    grad_psi: np.ndarray,
    *,
    lambda_v: float,
    v: float,
    a: float,
    b: float,
    c: float,
) -> np.ndarray:
    """Compute T_ab for the polynomial potential and kinetic functions.

    Gradients have trailing (t,r,phi,z) components. All other dimensions
    are broadcast sample dimensions. An arbitrary positive potential offset
    would change T_ab, so the potential is fixed to Eq. (14).
    """
    rho = np.asarray(rho, dtype=float)
    gradients = [np.asarray(g, dtype=float) for g in (grad_rho, grad_theta, grad_psi)]
    if any(g.shape[-1:] != (4,) for g in gradients):
        raise ValueError("each gradient must have four orthonormal components")
    if not np.isfinite(lambda_v) or lambda_v <= 0 or not np.isfinite(v) or v <= 0:
        raise ValueError("lambda_v and v must be finite and positive")
    if not all(np.all(np.isfinite(x)) for x in [rho, *gradients, a, b, c]):
        raise ValueError("all parameters and gradients must be finite")
    gr, gt, gp = gradients
    k_theta, k_psi, mixing = 1 + a * rho**2, 1 + b * rho**2, c * rho**2

    def dot(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.einsum("...i,ij,...j->...", x, ETA, y)

    potential = lambda_v / 4 * (rho**2 - v**2) ** 2
    lagrangian = (
        -0.5 * dot(gr, gr) - potential
        -0.5 * k_theta * dot(gt, gt) -0.5 * k_psi * dot(gp, gp)
        -mixing * dot(gt, gp)
    )

    def outer(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.einsum("...i,...j->...ij", x, y)

    return (
        outer(gr, gr)
        + k_theta[..., None, None] * outer(gt, gt)
        + k_psi[..., None, None] * outer(gp, gp)
        + mixing[..., None, None] * (outer(gt, gp) + outer(gp, gt))
        + lagrangian[..., None, None] * ETA
    )


def null_contraction(tensor: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """Contract T_ab with k=(1,n), requiring a unit spatial direction n."""
    direction = np.asarray(direction, dtype=float)
    if direction.shape[-1:] != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError("direction must have three finite components")
    if not np.allclose(np.linalg.norm(direction, axis=-1), 1.0, rtol=1e-10, atol=1e-10):
        raise ValueError("null direction must have unit spatial norm")
    k = np.concatenate((np.ones(direction.shape[:-1] + (1,)), direction), axis=-1)
    return np.einsum("...i,...ij,...j->...", k, tensor, k)
