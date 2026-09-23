"""Checks for the coupled matter evolution and its discrete conservation law."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from dee_radial_dynamics import initial_state, rhs, run, total_energy


def example_config():
    return {
        "outer_radius": 8.0, "width": 1.0, "grid_points": 120,
        "lambda_v": 4.0, "v": 1.0, "a": 0.3, "b": 0.2, "c": 0.1,
        "rho_amplitude": 0.02, "theta_amplitude": 0.2, "psi_amplitude": -0.1,
        "time_step": 0.0025, "steps": 400, "snapshot_interval": 200,
    }


class DynamicMatterTests(unittest.TestCase):
    def test_vacuum_is_stationary_with_zero_energy(self):
        config = example_config()
        for key in ("rho_amplitude", "theta_amplitude", "psi_amplitude"):
            config[key] = 0.0
        state = initial_state(config)
        np.testing.assert_array_equal(rhs(state, config), np.zeros_like(state))
        self.assertEqual(total_energy(state, config), 0.0)

    def test_rhs_preserves_discrete_energy_to_first_order(self):
        config = example_config()
        state = initial_state(config)
        r = (np.arange(config["grid_points"]) + 0.5) * config["outer_radius"] / config["grid_points"]
        state[3] = 0.03 * np.exp(-r**2)
        state[4] = 0.04 * np.exp(-r**2)
        state[5] = -0.02 * np.exp(-r**2)
        derivative = rhs(state, config)
        epsilon = 1e-5
        energy_rate = (total_energy(state + epsilon * derivative, config)
                       - total_energy(state - epsilon * derivative, config)) / (2 * epsilon)
        self.assertLess(abs(energy_rate), 1e-8)

    def test_nontrivial_run_conserves_energy_and_records_nec(self):
        config = example_config()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summary = run(config, output)
            self.assertEqual(len((output / "dynamic_profiles.csv").read_text().splitlines()), 361)
        self.assertLess(summary["max_relative_energy_drift"], 1e-8)
        self.assertLess(summary["snapshots"][-1]["minimum_rho"], 0.99)
        self.assertTrue(all(s["minimum_kinetic_determinant"] > 0 for s in summary["snapshots"]))
        self.assertTrue(all(s["minimum_radial_nec"] >= -1e-12 for s in summary["snapshots"]))
        self.assertTrue(all(s["minimum_azimuthal_nec"] >= -1e-12 for s in summary["snapshots"]))


if __name__ == "__main__":
    unittest.main()
