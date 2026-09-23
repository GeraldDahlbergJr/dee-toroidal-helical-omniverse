"""Checks of the action-derived local tensor and radial consistency gate."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from dee_matter import null_contraction, stress_energy
from matter_profile_check import run


class MatterTensorTests(unittest.TestCase):
    def test_null_contraction_matches_gradient_squares_in_mixed_case(self):
        rho = 1.1
        gr = np.array([0.0, 0.3, 0.0, 0.0])
        gt = np.array([0.4, 0.0, 0.2, 0.1])
        gp = np.array([-0.1, 0.05, -0.2, 0.0])
        a, b, c = 0.3, 0.2, 0.1
        tensor = stress_energy(rho, gr, gt, gp, lambda_v=4, v=1, a=a, b=b, c=c)
        for direction in (np.array([1.0, 0.0, 0.0]), np.array([0.0, -1.0, 0.0]),
                          np.array([0.3, 0.4, np.sqrt(0.75)])):
            k = np.r_[1.0, direction]
            u, w = np.dot(k, gt), np.dot(k, gp)
            expected = (np.dot(k, gr)**2 + (1 + a * rho**2) * u**2
                        + (1 + b * rho**2) * w**2 + 2 * c * rho**2 * u * w)
            self.assertAlmostEqual(null_contraction(tensor, direction), expected)

    def test_potential_contributes_only_metric_term(self):
        zero = np.zeros(4)
        tensor = stress_energy(0, zero, zero, zero, lambda_v=4, v=1, a=0, b=0, c=0)
        np.testing.assert_allclose(tensor, np.diag([1.0, -1.0, -1.0, -1.0]))
        self.assertAlmostEqual(null_contraction(tensor, np.array([0.0, 0.0, 1.0])), 0.0)

    def test_prescribed_profile_fails_conservation_gate(self):
        config = {
            "outer_radius": 6.0, "width": 1.0, "grid_points": 601,
            "lambda_v": 4.0, "v": 1.0, "a": 0.3, "b": 0.2, "c": 0.1,
            "rho_amplitude": 0.05, "theta_amplitude": 0.2, "psi_amplitude": -0.1,
        }
        with tempfile.TemporaryDirectory() as directory:
            summary = run(config, Path(directory))
            self.assertEqual(json.loads((Path(directory) / "summary.json").read_text()), summary)
        self.assertGreater(summary["max_abs_interior_residuals"]["conservation"], 0.01)
        self.assertGreater(summary["max_abs_interior_residuals"]["rho_equation"], 0.1)
        self.assertGreaterEqual(summary["minimum_radial_nec"], -1e-14)


if __name__ == "__main__":
    unittest.main()
