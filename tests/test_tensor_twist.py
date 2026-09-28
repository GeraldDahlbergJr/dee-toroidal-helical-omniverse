"""Independent analytic and numerical checks for the tensor twist sector."""

import unittest

import numpy as np

from tensor_twist import solve_tensor_twist


class TensorTwistTests(unittest.TestCase):
    def test_analytic_operator_identities(self):
        r = np.linspace(0.1, 2.0, 50)
        # L[H] = H'' - H'/r. Polynomial identities are exact analytically.
        for power, coefficient in ((2, 0.0), (4, 8.0), (6, 24.0)):
            expected = coefficient * r ** (power - 2)
            analytic = power * (power - 2) * r ** (power - 2)
            np.testing.assert_allclose(analytic, expected, rtol=0.0, atol=1e-14)

    def test_manufactured_solution_second_order_convergence(self):
        G = 0.2
        a = 0.3
        b = 0.7
        R = 2.0
        errors = []
        for points in (41, 81, 161, 321):
            r = np.linspace(0.0, R, points)
            exact = a * r**2 + b * r**4
            # L[exact] = 8 b r^2 = -16 pi G T_phi_z.
            stress = -(8.0 * b * r**2) / (16.0 * np.pi * G)
            numerical = solve_tensor_twist(r, stress, G, exact[-1])
            errors.append(np.max(np.abs(numerical - exact)))
            self.assertAlmostEqual(numerical[0], 0.0)
            self.assertAlmostEqual(numerical[-1], exact[-1])

        orders = [np.log(errors[i] / errors[i + 1]) / np.log(2.0) for i in range(3)]
        for order in orders:
            self.assertGreater(order, 1.9)

    def test_rejects_nonuniform_grid(self):
        r = np.array([0.0, 0.2, 0.5, 1.0])
        with self.assertRaisesRegex(ValueError, "uniform"):
            solve_tensor_twist(r, np.zeros_like(r))


if __name__ == "__main__":
    unittest.main()
