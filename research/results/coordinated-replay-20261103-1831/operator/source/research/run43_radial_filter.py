"""Original tensor constraint with a consistent sixth-difference radial filter.

The radial coefficient h^4/16 gives the same Nyquist damping -4/h^2
as the baseline angular fourth-difference filter. No field or gate fitting.
Eight-point D6 boundary closures are polynomial exact through degree seven.
Original sources, longitudinal tensor, scalar operator and audit are unchanged.
This retains the baseline O(h) first-interior-row composition defect; it is a
stability alternative, not a claim of uniform second-order boundary truncation.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator
import run43_consistent_momentum as baseline


def sixth(n,h):
    a=sp.lil_matrix((n,n))
    for i in range(n):
        start=min(max(i-3,0),n-8);js=np.arange(start,start+8)
        offsets=(js-i).astype(float);rhs=np.zeros(8);rhs[6]=720
        weights=np.linalg.solve(np.array([offsets**k for k in range(8)]),rhs)/h**6
        for j,w in zip(js,weights):a[i,j]+=w
    return a.tocsr()


class FilterGrid(baseline.TensorGrid):
    def __init__(self,n,stabilization=.25):
        super().__init__(n,stabilization)
        h=.2/n
        self.filter=sp.kron(sp.kron(sixth(n+1,h),sp.eye(n),format='csr'),sp.eye(n),format='csr')
        original=self.vector
        def apply(x):
            w=np.zeros((3,self.size));w[:,self.idx]=x.reshape(3,-1)
            v=(self.filter@w.T).T;er=self.basis[0]
            correction=h**4/16*(v+er*np.sum(er*v,axis=0)/3)
            return original@x+correction[:,self.idx].ravel()
        self.vector=LinearOperator(original.shape,matvec=apply,dtype=float)


if __name__=='__main__':
    solve=FilterGrid.solve
    def snapshot(self):
        out=baseline.output_dir()
        (out/'executed_radial_filter.py').write_bytes(Path(__file__).read_bytes())
        (out/'correction.json').write_text(json.dumps(dict(
            operator='original tensor divergence plus h^4/16 (I+er er^T/3) D6 W',
            source='unchanged',frozen_audit='unchanged',
            boundary='eight-point polynomial D6, no added physical boundary condition',
            limitation='retains original first-interior-row O(h) composition defect'),indent=2)+'\n')
        return solve(self)
    FilterGrid.solve=snapshot
    baseline.TensorGrid=FilterGrid
    baseline.main()
