import numpy as np

from nonlinear_initial_data import (
    flat_homogeneous_equilibrium_residual,
    hamiltonian_residual,
)


def test_vacuum_flat_slice_satisfies_hamiltonian_constraint():
    gamma = np.eye(3)
    K = np.zeros((3, 3))
    H = hamiltonian_residual(0.0, K, gamma, 0.0)
    assert abs(H) < 1e-14


def test_flat_slice_with_positive_dee_energy_fails_constraint():
    rho2 = 0.945
    rho = np.sqrt(rho2)
    H, E = flat_homogeneous_equilibrium_residual(rho)
    assert E > 0.0
    assert H < 0.0
    assert abs(H) > 1e-8


def test_constraint_is_satisfied_if_curvature_matches_energy_on_time_symmetric_slice():
    E = 0.125
    R3 = 16.0 * np.pi * E
    H = hamiltonian_residual(R3, np.zeros((3,3)), np.eye(3), E)
    assert abs(H) < 1e-14
