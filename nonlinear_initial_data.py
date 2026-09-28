"""3+1 initial-data diagnostics for the nonlinear DEE validation branch.

This module does not claim to evolve the full Einstein-matter system.  It
implements the ADM Hamiltonian and momentum constraints and the matter
projections needed to test candidate initial data derived from the archived
DEE action.
"""
from __future__ import annotations

import numpy as np


def potential(rho, lambda_v=4.0, v=1.0):
    return 0.25 * lambda_v * (rho * rho - v * v) ** 2


def kinetic_coefficients(rho, a=0.30, b=0.20, c=0.10):
    r2 = rho * rho
    return 1.0 + a*r2, 1.0 + b*r2, c*r2


def matter_energy_density(rho, pi_rho, grad_rho, pi_theta, grad_theta,
                          pi_psi, grad_psi, gamma_inv,
                          lambda_v=4.0, v=1.0, a=0.30, b=0.20, c=0.10):
    """Eulerian energy density E=n^mu n^nu T_munu for the DEE action."""
    kt, kp, lam = kinetic_coefficients(rho, a, b, c)
    def norm2(g): return np.einsum('...i,...ij,...j->...', g, gamma_inv, g)
    cross = np.einsum('...i,...ij,...j->...', grad_theta, gamma_inv, grad_psi)
    return (0.5*(pi_rho*pi_rho + norm2(grad_rho)) + potential(rho, lambda_v, v)
            + 0.5*kt*(pi_theta*pi_theta + norm2(grad_theta))
            + 0.5*kp*(pi_psi*pi_psi + norm2(grad_psi))
            + lam*(pi_theta*pi_psi + cross))


def matter_momentum_density(rho, pi_rho, grad_rho, pi_theta, grad_theta,
                            pi_psi, grad_psi, a=0.30, b=0.20, c=0.10):
    """S_i=-gamma_i^mu n^nu T_munu; sign follows Pi=-n^mu d_mu field."""
    kt, kp, lam = kinetic_coefficients(rho, a, b, c)
    return (pi_rho[..., None]*grad_rho
            + kt[..., None]*pi_theta[..., None]*grad_theta
            + kp[..., None]*pi_psi[..., None]*grad_psi
            + lam[..., None]*(pi_theta[..., None]*grad_psi
                              + pi_psi[..., None]*grad_theta))


def hamiltonian_residual(R3, Kij, gamma_inv, energy, G=1.0):
    """H = R3 + K^2-K_ij K^ij -16 pi G E."""
    K = np.einsum('...ij,...ij->...', gamma_inv, Kij)
    Kup = np.einsum('...ik,...jl,...kl->...ij', gamma_inv, gamma_inv, Kij)
    K2 = np.einsum('...ij,...ij->...', Kij, Kup)
    return R3 + K*K - K2 - 16.0*np.pi*G*energy


def flat_homogeneous_equilibrium_residual(rho_star, G=1.0, **matter):
    """Diagnostic only: flat, K_ij=0 data with nonzero DEE energy cannot satisfy H=0."""
    shape = np.shape(rho_star)
    gamma_inv = np.broadcast_to(np.eye(3), shape + (3,3))
    zeros_v = np.zeros(shape + (3,))
    zeros = np.zeros(shape)
    energy = matter_energy_density(rho_star, zeros, zeros_v, zeros, zeros_v,
                                   zeros, zeros_v, gamma_inv, **matter)
    Kij = np.zeros(shape + (3,3))
    return hamiltonian_residual(zeros, Kij, gamma_inv, energy, G), energy
