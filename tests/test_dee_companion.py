"""Automated checks for the DEE computational companion."""

import unittest
import warnings
import json
import tempfile
from pathlib import Path

import numpy as np

from dee_field_model import (
    effective_mass_squared,
    finite_equilibrium,
    finite_equilibrium_squared,
    is_positive_definite,
    kinetic_determinant,
    kinetic_matrix,
)
from linearized_einstein import (
    LinearizedCylindricalApproximationWarning,
    solve_cylindrical_linearized_field,
)
from weak_field_checkpoint import run
from nec_scan import null_energy_quadratic_form, scan_transport_projections


class DeeFieldModelTests(unittest.TestCase):
    def test_analytic_finite_equilibrium_and_mass(self):
        rho_squared = finite_equilibrium_squared(3.0, 2.0, 1.0, 0.5, 1.0, 2.0, 3.0, 4.0)
        self.assertAlmostEqual(rho_squared, 7.25)
        self.assertAlmostEqual(finite_equilibrium(3.0, 2.0, 1.0, 0.5, 1.0, 2.0, 3.0, 4.0), np.sqrt(7.25))
        self.assertAlmostEqual(effective_mass_squared(4.0, rho_squared), 58.0)

    def test_kinetic_matrix_positive_definiteness(self):
        rho_squared = 0.5
        matrix = kinetic_matrix(3.0, 2.0, 1.0, rho_squared)
        np.testing.assert_allclose(matrix, [[2.5, 0.5], [0.5, 2.0]])
        self.assertAlmostEqual(kinetic_determinant(3.0, 2.0, 1.0, rho_squared), 4.75)
        self.assertTrue(is_positive_definite(3.0, 2.0, 1.0, rho_squared))
        self.assertFalse(is_positive_definite(-2.0, -2.0, 0.0, 1.0))

    def test_rejects_invalid_rho_squared_inputs(self):
        for function, arguments in (
            (kinetic_matrix, (1.0, 1.0, 0.0, -0.1)),
            (kinetic_determinant, (1.0, 1.0, 0.0, np.nan)),
            (is_positive_definite, (1.0, 1.0, 0.0, np.inf)),
            (effective_mass_squared, (1.0, -0.1)),
        ):
            with self.subTest(function=function.__name__):
                with self.assertRaisesRegex(ValueError, "finite non-negative"):
                    function(*arguments)

    def test_rejects_invalid_lambda_v_inputs(self):
        for lambda_v in (0.0, -1.0, np.nan, np.inf):
            with self.subTest(function="finite_equilibrium_squared", lambda_v=lambda_v):
                with self.assertRaisesRegex(ValueError, "finite strictly positive"):
                    finite_equilibrium_squared(3.0, 2.0, 1.0, 0.5, 1.0, 2.0, 3.0, lambda_v)
            with self.subTest(function="effective_mass_squared", lambda_v=lambda_v):
                with self.assertRaisesRegex(ValueError, "finite strictly positive"):
                    effective_mass_squared(lambda_v, 1.0)


class NecTests(unittest.TestCase):
    def test_nec_quadratic_form_and_scan(self):
        matrix = kinetic_matrix(3.0, 2.0, 0.5, rho_squared=1.0)
        self.assertAlmostEqual(null_energy_quadratic_form(matrix, np.array([1.0, -2.0])), 14.0)
        scan = scan_transport_projections(matrix, samples=720)
        self.assertTrue(scan["all_nonnegative"])
        self.assertGreater(scan["minimum"], 0.0)

    def test_indefinite_matrix_has_negative_scan_direction(self):
        scan = scan_transport_projections(kinetic_matrix(0.0, -2.0, 0.0, rho_squared=1.0), samples=360)
        self.assertFalse(scan["all_nonnegative"])
        self.assertLess(scan["minimum"], 0.0)


class LinearizedEinsteinTests(unittest.TestCase):
    @staticmethod
    def _exact_field(radius, gravitational_constant=0.2, t0=1.3, alpha=0.7):
        source_scale = -16.0 * np.pi * gravitational_constant * t0
        outer = radius[-1]
        return source_scale * ((radius**2 - outer**2) / 4.0 + alpha * (radius**4 - outer**4) / 16.0)

    def test_radial_solver_grid_convergence(self):
        errors = []
        for points in (41, 81, 161):
            radius = np.linspace(0.0, 2.0, points)
            stress = 1.3 * (1.0 + 0.7 * radius**2)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", LinearizedCylindricalApproximationWarning)
                field = solve_cylindrical_linearized_field(radius, stress, gravitational_constant=0.2)
            errors.append(np.max(np.abs(field - self._exact_field(radius))))
        self.assertLess(errors[1], errors[0] / 3.5)
        self.assertLess(errors[2], errors[1] / 3.5)

    def test_regular_axis_and_outer_boundary(self):
        radius = np.linspace(0.0, 1.0, 101)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", LinearizedCylindricalApproximationWarning)
            field = solve_cylindrical_linearized_field(radius, np.ones_like(radius), outer_boundary=2.5)
        self.assertAlmostEqual(field[-1], 2.5)
        self.assertAlmostEqual(field[1] - field[0], -4.0 * np.pi * radius[1] ** 2, places=10)

    def test_rejects_nonfinite_radius_stress_and_boundary(self):
        radius = np.linspace(0.0, 1.0, 5)
        cases = (
            (np.array([0.0, 0.25, np.inf, 0.75, 1.0]), np.ones(5), 0.0, "radius"),
            (radius, np.array([1.0, 1.0, np.nan, 1.0, 1.0]), 0.0, "stress_energy"),
            (radius, np.ones(5), np.inf, "outer_boundary"),
        )
        for invalid_radius, invalid_stress, invalid_boundary, message in cases:
            with self.subTest(message=message), warnings.catch_warnings():
                warnings.simplefilter("ignore", LinearizedCylindricalApproximationWarning)
                with self.assertRaisesRegex(ValueError, message):
                    solve_cylindrical_linearized_field(
                        invalid_radius, invalid_stress, outer_boundary=invalid_boundary
                    )

    def test_rejects_incompatible_component_boundary_shape(self):
        radius = np.linspace(0.0, 1.0, 5)
        stress = np.ones((5, 2))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", LinearizedCylindricalApproximationWarning)
            with self.assertRaisesRegex(ValueError, "broadcastable"):
                solve_cylindrical_linearized_field(radius, stress, outer_boundary=np.ones(3))


class CheckpointTests(unittest.TestCase):
    def test_outputs_and_independent_exact_solution(self):
        config = {
            "outer_radius": 3.0,
            "source_amplitude": 0.01,
            "gravitational_constant": 1.0,
            "source_scale": 1e-6,
            "grid_points": [101, 201, 401],
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summary = run(config, output)
            self.assertEqual(json.loads((output / "summary.json").read_text()), summary)
            self.assertEqual(len((output / "radial_profile.csv").read_text().splitlines()), 402)
        errors = list(summary["max_absolute_error_by_grid"].values())
        self.assertLess(errors[2], errors[1] / 3.9)
        self.assertLess(errors[1], errors[0] / 3.9)
        self.assertLess(summary["source_scaling_absolute_difference"], 1e-13)


if __name__ == "__main__":
    unittest.main()
