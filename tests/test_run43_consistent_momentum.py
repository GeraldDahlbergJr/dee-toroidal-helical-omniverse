import sys, unittest, os
from unittest.mock import patch
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from run43_consistent_momentum import TensorGrid, scalar_flux_laplacian, FourierPreconditioner
from diagnose_toroidal_boundary_convergence import Audit

class ConsistentMomentum(unittest.TestCase):
    def setUp(self):
        # Production experiment settings must not select a different scalar
        # operator/preconditioner while testing the baseline algebra.
        self.environment = patch.dict(os.environ)
        self.environment.start()
        for key in ('DEE_SCALAR_ORDER', 'DEE_PRECONDITIONER'):
            os.environ.pop(key, None)
        self.addCleanup(self.environment.stop)

    def test_fourier_inverse_and_scalar_flux_consistency(self):
        g=TensorGrid(8);L=scalar_flux_laplacian(g);interior=(L[g.idx]@g.inject).tocsc()
        pre=FourierPreconditioner(interior,8)
        rng=np.random.default_rng(43);rhs=rng.normal(size=len(g.idx))
        np.testing.assert_allclose(interior@pre.solve(rhs),rhs,atol=1e-12,rtol=1e-11)
        u=1+1e-3*rng.normal(size=g.size);w=np.zeros((3,)+g.shape)
        audit=Audit(g.shape);h=audit.flux_hamiltonian(u.reshape(g.shape),w).ravel()
        # Remove exactly the declared angular stabilization and compare fluxes.
        raw=u.reshape(g.shape);ht=2*np.pi/8
        stab=np.zeros(g.size)
        for axis,scale in [(1,g.r**2),(2,g.q**2)]:
            d4=np.roll(raw,2,axis)-4*np.roll(raw,1,axis)+6*raw-4*np.roll(raw,-1,axis)+np.roll(raw,-2,axis)
            stab+=.25*d4.ravel()/ht**2/scale
        expected=-8*u**-5*(L@u+stab)-16*np.pi*1e-3*g.energy(u)
        np.testing.assert_allclose(expected[g.idx],h[g.idx],atol=1e-10,rtol=1e-11)

    def test_tensor_operator_matches_independent_audit(self):
        g=TensorGrid(8,stabilization=0)
        rng=np.random.default_rng(43);x=rng.normal(size=g.vector.shape[0])
        w=np.zeros((3,g.size));w[:,g.idx]=x.reshape(3,-1)
        a=Audit(g.shape);tensor=a.tensor(w.reshape((3,)+g.shape),2)
        expected=np.array([sum(a.grad(tensor[i,j],2)[j] for j in range(3)) for i in range(3)]).reshape(3,-1)[:,g.idx].ravel()
        np.testing.assert_allclose(g.vector@x,expected,rtol=1e-12,atol=1e-9)

    def test_stabilization_damps_angular_checkerboard(self):
        g=TensorGrid(8);raw=TensorGrid(8,stabilization=0)
        radial=np.sin(np.pi*(g.r-.05)/.20)
        pattern=(-1.)**np.indices(g.shape)[1].ravel()
        w=np.array([radial*pattern,np.zeros(g.size),np.zeros(g.size)])
        x=w[:,g.idx].ravel();delta=(g.vector@x)-(raw.vector@x)
        self.assertLess(float(x@delta),0)
        self.assertGreater(np.linalg.norm(g.vector@x)/np.linalg.norm(x),1)

    def test_small_grid_has_no_exact_null_mode(self):
        # Full singular spectrum also detects mixed and radial checkerboards.
        g=TensorGrid(6);eye=np.eye(g.vector.shape[0])
        matrix=np.column_stack([g.vector@x for x in eye])
        s=np.linalg.svd(matrix,compute_uv=False)
        self.assertGreater(s[-1]/s[0],1e-7)

if __name__=='__main__':unittest.main()
