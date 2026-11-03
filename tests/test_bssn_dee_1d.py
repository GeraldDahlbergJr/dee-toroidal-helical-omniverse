"""Analytic conformal-curvature and coupled-evolution checks."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
from evolve_dee_adm_1d import L, geometry, initial as adm_initial
from evolve_dee_bssn_variables_1d import (conformal_ricci,to_bssn,to_adm,run)


class BSSNTest(unittest.TestCase):
    def test_conformal_ricci_matches_analytic_conformal_metric(self):
        n=256; h=L/n; x=np.arange(n)*h; omega=2*np.pi/L
        phi=.02*np.sin(omega*x)
        tg=np.ones((n,3)); conn=np.zeros(n)
        ric=conformal_ricci(phi,tg,conn,h)
        second=-.02*omega**2*np.sin(omega*x)
        first=.02*omega*np.cos(omega*x)
        expected=np.column_stack((-4*second,-2*second-4*first**2,
                                  -2*second-4*first**2))
        np.testing.assert_allclose(ric,expected,atol=3e-6)

    def test_evolved_connection_changes_curvature(self):
        n=128; h=L/n; x=np.arange(n)*h
        tg=np.ones((n,3)); phi=np.zeros(n)
        perturb=.01*np.sin(2*np.pi*x/L)
        ric=conformal_ricci(phi,tg,perturb,h)
        np.testing.assert_allclose(ric[:,0],np.gradient(perturb,h,edge_order=2),atol=3e-4)
        np.testing.assert_allclose(ric[:,1:],0,atol=1e-14)

    def test_round_trip_and_short_evolution(self):
        _,adm=adm_initial(64)
        for left,right in zip(adm,to_adm(to_bssn(adm,L/64))):
            np.testing.assert_allclose(left,right,rtol=1e-13,atol=1e-13)
        state,row=run(64)
        self.assertGreater(row['final']['metric_min'],0)
        self.assertLess(row['final']['conformal_det_max'],1e-12)
        self.assertLess(row['final']['conformal_trace_max'],1e-12)
        self.assertLess(row['final']['connection_constraint_l2'],1e-7)


if __name__=='__main__': unittest.main()
