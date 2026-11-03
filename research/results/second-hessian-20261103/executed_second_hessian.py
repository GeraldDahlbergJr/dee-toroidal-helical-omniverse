"""Standard second-order direct chart Hessian with original physical sources.

Direct radial D2 at every interior row removes composed endpoint error.
All angular, mixed and connection terms have the same second-order accuracy.
No vector dissipation, source correction, or fitted parameter is used.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
import run43_fourth_vector as fourth
import run43_consistent_momentum as baseline


def second_stencil(n,h,degree,periodic=False):
    count=3 if degree==1 or periodic else 4
    a=sp.lil_matrix((n,n))
    for i in range(n):
        start=min(max(i-1,0),n-count)
        js=np.arange(i-1,i+2) if periodic else np.arange(start,start+count)
        offsets=(js-i).astype(float);rhs=np.zeros(count);rhs[degree]=1 if degree==1 else 2
        weights=np.linalg.solve(np.array([offsets**k for k in range(count)]),rhs)/h**degree
        for j,w in zip(js,weights):a[i,j%n]+=w
    return a.tocsr()


class SecondHessianGrid(fourth.FourthVectorGrid):
    def __init__(self,n,stabilization=.25):
        # Reuse the analytic tensor assembly, substituting only derivative
        # matrices. This process-local assignment is restored before returning.
        original=fourth.stencil
        fourth.stencil=second_stencil
        try:super().__init__(n,stabilization)
        finally:fourth.stencil=original


if __name__=='__main__':
    solve=SecondHessianGrid.solve
    def snapshot(self):
        out=baseline.output_dir()
        (out/'executed_second_hessian.py').write_bytes(Path(__file__).read_bytes())
        (out/'executed_fourth_vector.py').write_bytes(Path(fourth.__file__).read_bytes())
        (out/'effective-settings.json').write_text(json.dumps(dict(
            vector_order=2,vector_stabilization=0,scalar_order=2,scalar_stabilization=.25,
            radial_D1_points=3,radial_D2_points=4,periodic_D1_D2_points=3,
            source='unchanged original',frozen_audit='unchanged'),indent=2)+'\n')
        return solve(self)
    SecondHessianGrid.solve=snapshot
    baseline.TensorGrid=SecondHessianGrid
    baseline.main()
