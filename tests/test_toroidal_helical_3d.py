"""Independent geometry and regularity checks for the 3-D trial."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
from toroidal_helical_initial_data import fields,solve,lap,L
from evolve_dee_ring_3d import geometric,run


class RingTrialTest(unittest.TestCase):
    def test_conformal_initial_constraint_and_smooth_axis(self):
        f,psi,k,record=solve(24)
        self.assertLess(record['hamiltonian_max'],1e-10)
        self.assertEqual(record['momentum_max'],0)
        # The integer angular mode is represented by a Cartesian polynomial:
        # no branch-cut angular coordinate or singular field on the axis.
        self.assertTrue(np.all(np.isfinite(f)))
        np.testing.assert_allclose(f[12,12,:,1:],0,atol=1e-15)

    def test_ricci_conformal_identity_refines(self):
        errors=[]
        for n in (12,24):
            h=L/n;x=np.arange(n)*h
            psi=1+.01*np.sin(2*np.pi*x[:,None,None]/L)*np.ones((1,n,n))
            gamma=np.zeros((n,n,n,3,3));curvature=np.zeros_like(gamma)
            for i in range(3):gamma[...,i,i]=psi**4
            _,_,scalar,_,_=geometric(gamma,curvature,h)
            exact=-8*psi**-5*lap(psi,h)
            errors.append(float(np.sqrt(np.mean((scalar-exact)**2))))
        self.assertGreater(errors[0]/errors[1],3.5)

    def test_off_diagonal_geometry_generated(self):
        state,record=run(12)
        self.assertGreater(np.max(np.abs(state[1][...,0,1])),1e-7)
        self.assertGreater(np.max(np.abs(state[0][...,0,1])),1e-8)
        self.assertGreater(record['final']['metric_eigenvalue_min'],0)

    def test_momentum_operator_against_flat_manufactured_data(self):
        n=24;h=L/n;x=np.arange(n)*h
        metric=np.broadcast_to(np.eye(3),(n,n,n,3,3)).copy()
        k=np.zeros_like(metric)
        shape=np.sin(2*np.pi*x/L)[:,None,None]*np.ones((1,n,n))
        k[...,0,0]=shape
        k[...,1,1]=np.cos(2*np.pi*x/L)[:,None,None]
        _,_,_,momentum,_=geometric(metric,k,h)
        expected=(np.sin(2*np.pi*h/L)/h)*np.sin(2*np.pi*x/L)[:,None,None]*np.ones((1,n,n))
        np.testing.assert_allclose(momentum[...,0],expected,atol=1e-13)
        np.testing.assert_allclose(momentum[...,1:],0,atol=1e-13)


if __name__=='__main__': unittest.main()
