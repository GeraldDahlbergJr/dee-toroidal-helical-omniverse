"""Radial principal-part consistency correction, without endpoint retuning.

Replace the composed radial second derivative by the direct three-point Drr
in the principal tensor term (I+er er^T/3). This is the radial counterpart of
the baseline's angular checkerboard correction. No adjustable coefficient.
"""
import sys,json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator
import run43_consistent_momentum as baseline
from solve_toroidal_helical_coupled_constraints import derivative
from validate_toroidal_helical_dee_source import RMIN,RMAX


class RadialGrid(baseline.TensorGrid):
    def __init__(self,n,stabilization=.25):
        super().__init__(n,stabilization)
        h=(RMAX-RMIN)/n
        drr=sp.kron(sp.kron(derivative(n+1,h,second=True),sp.eye(n),format='csr'),sp.eye(n),format='csr')
        self.radial_defect=(drr-self.d[0]@self.d[0]).tocsr()
        original=self.vector
        def apply(x):
            full=np.zeros((3,self.size));full[:,self.idx]=x.reshape(3,-1)
            defect=np.array([self.radial_defect@v for v in full])
            er=self.basis[0]
            correction=defect+er*np.sum(er*defect,axis=0)/3
            return original@x+correction[:,self.idx].ravel()
        self.vector=LinearOperator(original.shape,matvec=apply,dtype=float)

    def solve_vector(self,u):
        w=super().solve_vector(u)
        rhs=(8*np.pi*baseline.G*u[self.idx]**6*self.S[:,self.idx]).ravel()
        residual=float(np.linalg.norm(self.vector@w[:,self.idx].ravel()-rhs))
        tol=max(1e-13,1e-10*float(np.linalg.norm(rhs)))
        if not np.isfinite(residual) or residual>tol:
            raise RuntimeError(f'Vector residual {residual} exceeds unchanged tolerance {tol}')
        self.history[-1]['acceptance_tolerance']=tol
        return w


if __name__=='__main__':
    original_solve=RadialGrid.solve
    def solve_with_snapshot(self):
        out=baseline.output_dir()
        (out/'executed_radial_correction.py').write_bytes(Path(__file__).read_bytes())
        (out/'correction.json').write_text(json.dumps(dict(
            correction='(I+er er^T/3)(Drr-Dr Dr)W',
            coefficient='fixed by continuum principal symbol, no free parameter',
            physical_source='unchanged',frozen_audit='unchanged'),indent=2)+'\n')
        return original_solve(self)
    RadialGrid.solve=solve_with_snapshot
    baseline.TensorGrid=RadialGrid
    baseline.main()
