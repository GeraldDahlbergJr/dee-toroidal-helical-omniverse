"""Independent analytic manufactured checks for toroidal differential operators."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
from solve_toroidal_helical_coupled_constraints import Grid
from dee_stress_energy_3p1 import projections
from validate_toroidal_helical_dee_source import RHO,PI,KTH,KPS,M,metric

class ToroidalOperators(unittest.TestCase):
    def test_cartesian_polynomial_vector_laplacian_converges(self):
        errors=[]
        for n in (8,16):
            g=Grid(n+1,n,n)
            x=g.q*np.cos(g.p); y=g.q*np.sin(g.p); z=g.r*np.sin(g.t)
            w=np.array([x*x+y*y+z*z,x*y,z*z])
            div=sum(g.grad[i]@w[i] for i in range(3))
            actual=np.array([g.lap@w[i]+g.grad[i]@div/3 for i in range(3)])
            truth=np.array([7.,0.,8/3])[:,None]
            bulk=(g.r>=.1)&(g.r<=.2)
            errors.append(np.sqrt(np.mean((actual[:,bulk]-truth)**2)))
        self.assertGreater(np.log2(errors[0]/errors[1]),1.5)

    def test_source_uses_physical_inverse_metric(self):
        g=Grid(9,8,8); u=np.full(g.size,1.07)
        for n in (0,25,100):
            chart=metric(g.r[n],g.t[n]); D=np.zeros((3,3))
            D[1]=KTH*np.array([0.,1.,-float(M)])
            D[2]=KPS*np.array([0.,1.,-float(M)])
            E,S,_,_=projections(RHO,PI,D,np.linalg.inv(chart)/u[n]**4)
            self.assertAlmostEqual(E,g.energy(u)[n],places=12)
            er,et,ep=g.basis
            cart=er[:,n]*S[0]+et[:,n]*S[1]/g.r[n]+ep[:,n]*S[2]/g.q[n]
            np.testing.assert_allclose(cart,g.S[:,n],atol=1e-13)

if __name__=='__main__': unittest.main()
