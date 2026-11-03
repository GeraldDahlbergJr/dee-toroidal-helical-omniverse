"""Matched radial endpoint experiment; original matter and frozen audit unchanged.

Match the +h^2 f'''/6 error of centered D at both endpoints. This avoids
the O(h) error in the first interior row of D(Df). No gate parameter is tuned.
The unchanged frozen audit continues to use its original endpoint derivatives.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
import run43_consistent_momentum as baseline
from investigate_toroidal_inner_operator import endpoint_cubic_matrix
from validate_toroidal_helical_dee_source import RMIN, RMAX


class EndpointGrid(baseline.TensorGrid):
    def __init__(self, n, stabilization=.25):
        super().__init__(n, stabilization)
        dr = sp.csr_matrix(endpoint_cubic_matrix(n+1, (RMAX-RMIN)/n, matched=True))
        self.d[0] = sp.kron(sp.kron(dr, sp.eye(n), format='csr'), sp.eye(n), format='csr')
        er, et, ep = self.basis
        self.grad = [(sp.diags(er[i]) @ self.d[0]
                      + sp.diags(et[i]/self.r) @ self.d[1]
                      + sp.diags(ep[i]/self.q) @ self.d[2]).tocsr() for i in range(3)]

    def solve_vector(self, u):
        w = super().solve_vector(u)
        rhs = (8*np.pi*baseline.G*u[self.idx]**6*self.S[:,self.idx]).ravel()
        tolerance = max(1e-13, 1e-10*np.linalg.norm(rhs))
        residual = np.linalg.norm(self.vector @ w[:,self.idx].ravel()-rhs)
        # Check the actual unpreconditioned residual, including nonfinite values.
        if not np.isfinite(residual) or residual > tolerance:
            raise RuntimeError(f'Vector residual {residual} exceeds {tolerance}')
        self.history[-1]['acceptance_tolerance'] = float(tolerance)
        return w


if __name__ == '__main__':
    # Reuse the exact source, nonlinear equations, tolerances, fields and frozen
    # twelve-check reporting. The additional snapshot identifies this correction.
    original_solve = EndpointGrid.solve
    def solve_with_snapshot(self):
        out = baseline.output_dir()
        (out/'executed_endpoint_correction.py').write_bytes(Path(__file__).read_bytes())
        (out/'correction.json').write_text(json.dumps(dict(
            radial='matched leading centered first-derivative truncation at endpoints',
            source='unchanged original Einstein-DEE projections',
            frozen_audit='unchanged', stabilization=.25), indent=2)+'\n')
        return original_solve(self)
    EndpointGrid.solve = solve_with_snapshot
    baseline.TensorGrid = EndpointGrid
    baseline.main()
