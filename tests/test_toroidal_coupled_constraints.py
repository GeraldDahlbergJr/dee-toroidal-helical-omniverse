"""Independent analytic manufactured checks for toroidal differential operators."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
from solve_toroidal_helical_coupled_constraints import Grid
from diagnose_toroidal_boundary_convergence import Audit, radial_matrix
from investigate_toroidal_inner_operator import endpoint_cubic_matrix, expanded
from solve_toroidal_matched_endpoint_candidate import CandidateGrid
from solve_toroidal_direct_hessian_candidate import DirectHessianGrid, manufactured_check
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

    def test_fourth_order_radial_polynomial_including_endpoints(self):
        x=np.linspace(.05,.25,9)
        derivative=radial_matrix(9,x[1]-x[0],4)
        for degree in range(5):
            expected=np.zeros_like(x) if degree==0 else degree*x**(degree-1)
            np.testing.assert_allclose(derivative@x**degree,expected,atol=2e-13)

    def test_independent_tensor_divergence_manufactured_vector(self):
        errors=[]
        for n in (8,16):
            a=Audit((n+1,n,n))
            x=a.q*np.cos(a.p); y=a.q*np.sin(a.p)
            z=np.broadcast_to(a.r*np.sin(a.t),a.shape)
            w=np.array([x*x+y*y+z*z,x*y,z*z])
            tensor=a.tensor(w,4)
            actual=np.array([sum(a.grad(tensor[i,j],4)[j] for j in range(3)) for i in range(3)])
            truth=np.array([7.,0.,8/3])[:,None,None,None]
            errors.append(np.sqrt(np.mean((actual-truth)**2)))
        self.assertGreater(np.log2(errors[0]/errors[1]),3.)

    def test_matched_endpoint_composition_has_no_cubic_edge_defect(self):
        x=np.linspace(.05,.25,17);h=x[1]-x[0]
        matched=endpoint_cubic_matrix(len(x),h,matched=True)
        np.testing.assert_allclose((matched@matched@x**3)[1:-1],6*x[1:-1],atol=2e-12)
        standard=radial_matrix(len(x),h,2)
        self.assertGreater(abs((standard@standard@x**3)[1]-6*x[1]),h/2)

    def test_expanded_audit_matches_production_vector_operator(self):
        g=Grid(9,8,8);a=Audit(g.shape)
        f=(g.r-.05)*(.25-g.r)*np.exp(8*(g.r-.05))
        w=np.array([.3*f,-.4*f,.5*f])
        expected=(g.vector@w[:,g.idx].ravel()).reshape(3,-1)
        actual=expanded(a,w.reshape((3,)+g.shape)).reshape(3,-1)[:,g.idx]
        np.testing.assert_allclose(actual,expected,atol=2e-11,rtol=1e-10)

    def test_candidate_preserves_source_and_matches_independent_operator(self):
        base=Grid(9,8,8);g=CandidateGrid(9,8,8);a=Audit(g.shape)
        a.drad[2]=endpoint_cubic_matrix(9,a.dr,matched=True)
        np.testing.assert_array_equal(g.S,base.S)
        np.testing.assert_array_equal(g.energy(np.ones(g.size)),base.energy(np.ones(g.size)))
        np.testing.assert_allclose((g.L-base.L).toarray(),0,atol=1e-12)
        f=(g.r-.05)*(.25-g.r)*np.exp(8*(g.r-.05))
        w=np.array([.3*f,-.4*f,.5*f])
        actual=expanded(a,w.reshape((3,)+g.shape)).reshape(3,-1)[:,g.idx]
        expected=(g.vector@w[:,g.idx].ravel()).reshape(3,-1)
        np.testing.assert_allclose(actual,expected,atol=2e-11,rtol=1e-10)

    def test_direct_hessian_manufactured_inner_boundary_converges(self):
        result=manufactured_check()
        self.assertTrue(result['gate'],result['orders'])

    def test_direct_hessian_polynomial_converges_and_preserves_matter(self):
        errors=[]
        for n in (8,16):
            g=DirectHessianGrid(n+1,n,n)
            x=g.q*np.cos(g.p);y=g.q*np.sin(g.p);z=g.r*np.sin(g.t)
            f=x*x+y*y+z*z
            actual=sum(sum(g._cartesian_hessian_coefficients[i][i][k]*(g._hessian_components[k]@f) for k in range(6)) for i in range(3))
            errors.append(np.sqrt(np.mean((actual[g.idx]-6.)**2)))
            a=Audit(g.shape)
            np.testing.assert_allclose(g.S,a.S.reshape(3,-1),atol=1e-14)
            np.testing.assert_allclose(g.energy(np.ones(g.size)),(.5*a.temporal+.5*a.spatial+g.U).ravel(),atol=1e-14)
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
