"""Validation gates for the action-derived nonlinear DEE matter sector.

These tests do NOT claim a nonlinear Einstein evolution has been solved. They
check exact algebraic consequences that must hold before constructing
constraint-satisfying 3+1 initial data.
"""

import unittest
import numpy as np

from dee_field_model import (
    effective_mass_squared,
    finite_equilibrium_squared,
    is_positive_definite,
    kinetic_matrix,
)


class NonlinearDEEMatterGate(unittest.TestCase):
    def setUp(self):
        self.lam = 4.0
        self.v = 1.0
        self.a = 0.30
        self.b = 0.20
        self.c = 0.10
        self.x_theta = 0.60
        self.x_psi = 0.40
        self.x_cross = -0.20
        self.rho2 = finite_equilibrium_squared(
            self.v, self.a, self.b, self.c,
            self.x_theta, self.x_psi, self.x_cross, self.lam
        )

    def test_archived_illustrative_point(self):
        self.assertAlmostEqual(self.rho2, 0.945, places=12)
        self.assertAlmostEqual(np.sqrt(self.rho2), np.sqrt(0.945), places=12)

    def test_finite_branch_exact_stationarity_residual(self):
        residual = self.lam * (self.rho2 - self.v**2)
        residual += self.a * self.x_theta
        residual += self.b * self.x_psi
        residual += 2.0 * self.c * self.x_cross
        self.assertAlmostEqual(residual, 0.0, places=13)

    def test_positive_radial_curvature_on_positive_branch(self):
        self.assertGreater(self.rho2, 0.0)
        self.assertGreater(effective_mass_squared(self.lam, self.rho2), 0.0)
        self.assertAlmostEqual(
            effective_mass_squared(self.lam, self.rho2),
            2.0 * self.lam * self.rho2,
            places=13,
        )

    def test_kinetic_matrix_is_positive_definite_at_archived_point(self):
        self.assertTrue(is_positive_definite(self.a, self.b, self.c, self.rho2))
        eig = np.linalg.eigvalsh(kinetic_matrix(self.a, self.b, self.c, self.rho2))
        self.assertTrue(np.all(eig > 0.0))

    def test_nec_quadratic_form_all_random_projected_directions(self):
        # The analytic result is eigenvalue positivity; random sampling is an
        # independent numerical regression guard, not the proof itself.
        K = kinetic_matrix(self.a, self.b, self.c, self.rho2)
        rng = np.random.default_rng(20260927)
        u = rng.normal(size=(100000, 2))
        q = np.einsum("ni,ij,nj->n", u, K, u)
        self.assertGreaterEqual(float(q.min()), -1e-13)

    def test_indefinite_kinetic_matrix_is_rejected(self):
        # Deliberately choose mixing large enough to make det(K)<0.
        rho2 = 1.0
        a, b, c = 0.1, 0.1, 2.0
        self.assertFalse(is_positive_definite(a, b, c, rho2))
        eig = np.linalg.eigvalsh(kinetic_matrix(a, b, c, rho2))
        self.assertLess(float(eig.min()), 0.0)


if __name__ == "__main__":
    unittest.main()
