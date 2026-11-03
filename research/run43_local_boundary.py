"""Local removal of the first-order radial composition defect.

Starting from the preserved tensor baseline, use direct Drr only at the two
adjacent interior planes. The continuum radial principal symbol fixes the
coefficient I+er er^T/3. No source, boundary datum or audit is changed.
This is a second-order local consistency repair, not a fourth-order claim.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator
import run43_consistent_momentum as baseline
from solve_toroidal_helical_coupled_constraints import derivative
from validate_toroidal_helical_dee_source import RMIN, RMAX


class LocalBoundaryGrid(baseline.TensorGrid):
    def __init__(self, n, stabilization=.25):
        super().__init__(n, stabilization)
        h = (RMAX-RMIN)/n
        drr = sp.kron(sp.kron(derivative(n+1, h, second=True),
                             sp.eye(n), format='csr'), sp.eye(n), format='csr')
        rows = np.zeros(self.shape)
        rows[[1, -2]] = 1
        self.radial_defect = (sp.diags(rows.ravel()) @
                             (drr-self.d[0]@self.d[0])).tocsr()
        original = self.vector
        def apply(x):
            full = np.zeros((3, self.size))
            full[:, self.idx] = x.reshape(3, -1)
            defect = (self.radial_defect @ full.T).T
            er = self.basis[0]
            correction = defect + er*np.sum(er*defect, axis=0)/3
            return original@x + correction[:, self.idx].ravel()
        self.vector = LinearOperator(original.shape, matvec=apply, dtype=float)

    def solve_vector(self, u):
        w = super().solve_vector(u)
        rhs = (8*np.pi*baseline.G*u[self.idx]**6*self.S[:, self.idx]).ravel()
        residual = float(np.linalg.norm(self.vector@w[:, self.idx].ravel()-rhs))
        tol = max(1e-13, 1e-10*float(np.linalg.norm(rhs)))
        if not np.isfinite(residual) or residual > tol:
            raise RuntimeError(f'True vector residual {residual} exceeds {tol}')
        self.history[-1]['acceptance_tolerance'] = tol
        return w


if __name__ == '__main__':
    solve = LocalBoundaryGrid.solve
    def snapshot(self):
        out = baseline.output_dir()
        (out/'executed_local_boundary.py').write_bytes(Path(__file__).read_bytes())
        (out/'correction.json').write_text(json.dumps(dict(
            correction='P_edge (I+er er^T/3)(Drr-Dr Dr)W',
            planes='first and last interior radial planes',
            coefficient='continuum principal symbol; fixed before execution',
            source='original unchanged', audit='frozen unchanged',
            extra_boundary_data=False, vector_boundary_consistency_order=2,
            angular_stabilization=.25, scalar='unchanged stabilized flux'),
            indent=2)+'\n')
        return solve(self)
    LocalBoundaryGrid.solve = snapshot
    baseline.TensorGrid = LocalBoundaryGrid
    baseline.main()
