import os
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'research'))
from run43_local_boundary import LocalBoundaryGrid
from run43_consistent_momentum import TensorGrid


class LocalBoundaryTests(unittest.TestCase):
    def test_cubic_defect_coefficient_and_locality(self):
        with patch.dict(os.environ, {'DEE_SCALAR_ORDER': 'flux',
                                    'DEE_PRECONDITIONER': 'fourier'}):
            g = LocalBoundaryGrid(8)
            b = TensorGrid(8)
        h = .2/8
        # Analytic cubic: direct Drr is exact; composition defect is +/-3h/2.
        v = g.r**3
        correction = g.radial_defect@v
        expected = np.zeros(g.shape)
        expected[1] = -1.5*h
        expected[-2] = 1.5*h
        np.testing.assert_allclose(correction.reshape(g.shape), expected,
                                   atol=2e-13, rtol=1e-10)
        rng = np.random.default_rng(43)
        x = rng.normal(size=g.vector.shape[0])
        delta = (g.vector@x-b.vector@x).reshape(3, 7, 8, 8)
        np.testing.assert_allclose(delta[:, 1:-1], 0, atol=1e-10)
        np.testing.assert_array_equal(g.S, b.S)
        np.testing.assert_array_equal(g.lap.toarray(), b.lap.toarray())
        u = np.linspace(1, 1.01, g.size)
        np.testing.assert_array_equal(g.energy(u), b.energy(u))


if __name__ == '__main__':
    unittest.main()
