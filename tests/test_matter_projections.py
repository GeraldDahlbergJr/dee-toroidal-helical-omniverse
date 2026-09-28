"""Independent covariant-tensor check of the DEE Eulerian projections."""
import unittest

import numpy as np

from nonlinear_initial_data import matter_projections


class ProjectionTest(unittest.TestCase):
    def test_mixed_kinetic_with_nontrivial_lapse_shift_metric(self):
        rho = 0.83
        lapse = 1.31
        shift = np.array([0.17, -0.12, 0.06])
        gamma = np.array([[1.7, 0.2, 0.1], [0.2, 1.3, -0.08],
                          [0.1, -0.08, 1.1]])
        inv = np.linalg.inv(gamma)
        pis = np.array([0.21, -0.34, 0.29])
        grads = np.array([[0.2, 0.1, -0.3], [0.4, -0.2, 0.1],
                          [-0.1, 0.5, 0.3]])
        a, b, c = 0.3, 0.2, 0.37
        e, s, sij, trace = matter_projections(
            rho, pis[0], grads[0], pis[1], grads[1], pis[2], grads[2],
            gamma, inv, a=a, b=b, c=c)
        # Construct the four-metric, coordinate gradients, and T_mu nu
        # directly from the covariant action, independently of the 3+1 formulas.
        g = np.zeros((4, 4))
        g[0, 0] = -lapse*lapse + shift @ gamma @ shift
        g[0, 1:] = gamma @ shift
        g[1:, 0] = g[0, 1:]
        g[1:, 1:] = gamma
        d = np.column_stack((-lapse*pis + grads @ shift, grads))
        field_metric = np.array([[1., 0., 0.],
                                 [0., 1+a*rho*rho, c*rho*rho],
                                 [0., c*rho*rho, 1+b*rho*rho]])
        contraction = np.einsum('ab,am,mn,bn->', field_metric, d,
                                np.linalg.inv(g), d)
        lagrangian = -0.5*contraction - (rho*rho-1)**2
        t = np.einsum('ab,am,bn->mn', field_metric, d, d) + g*lagrangian
        normal = np.r_[1., -shift]/lapse
        np.testing.assert_allclose(e, normal @ t @ normal, atol=1e-13)
        np.testing.assert_allclose(s, -normal @ t[:, 1:], atol=1e-13)
        np.testing.assert_allclose(sij, t[1:, 1:], atol=1e-13)
        np.testing.assert_allclose(trace, np.sum(inv*sij), atol=1e-13)

    def test_isotropic_potential(self):
        gamma = np.eye(3)
        z = np.zeros(3)
        e, s, sij, trace = matter_projections(
            0., 0., z, 0., z, 0., z, gamma, gamma)
        self.assertAlmostEqual(e, 1.)
        np.testing.assert_allclose(s, z)
        np.testing.assert_allclose(sij, -gamma)
        self.assertAlmostEqual(trace, -3.)


if __name__ == '__main__':
    unittest.main()
