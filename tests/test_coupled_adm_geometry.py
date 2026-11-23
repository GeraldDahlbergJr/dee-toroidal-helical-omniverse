"""Analytic geometry and constraint checks for the restricted ADM solver."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'reproduce'))
from evolve_dee_adm_1d import geometry, initial, diagnostics, L


class GeometryTest(unittest.TestCase):
    def test_conformal_scalar_curvature_identity(self):
        n=512; h=L/n; x=np.arange(n)*h
        psi=1+.03*np.sin(2*np.pi*x/L)
        gamma=np.repeat(psi[:,None]**4,3,axis=1)
        ric, ham, mom, _=geometry(gamma,np.zeros_like(gamma),h)
        exact=-8*psi**-5*(-.03*(2*np.pi/L)**2*np.sin(2*np.pi*x/L))
        self.assertLess(np.max(np.abs(ham-exact)), 2e-5)
        self.assertLess(np.max(np.abs(mom)), 1e-13)
        self.assertTrue(np.all(np.isfinite(ric)))

    def test_constraint_initial_data_refines(self):
        residuals=[]
        for n in (64,128,256):
            _,state=initial(n)
            residuals.append(diagnostics(state,L/n)['hamiltonian_l2'])
        self.assertGreater(np.log2(residuals[0]/residuals[1]),1.9)
        self.assertGreater(np.log2(residuals[1]/residuals[2]),1.9)


if __name__ == '__main__': unittest.main()
