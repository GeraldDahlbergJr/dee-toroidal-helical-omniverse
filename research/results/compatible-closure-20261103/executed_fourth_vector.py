"""Fourth-order direct Hessian momentum candidate, fixed before gate evaluation.

Five-point first and periodic second derivatives; six-point radial D2 closures.
Mixed chart derivatives and analytic toroidal connection terms construct the
Cartesian Hessian directly. Original matter, scalar flux and audit unchanged.
No vector stabilization is needed for the direct elliptic second derivatives.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator
import run43_consistent_momentum as baseline
from validate_toroidal_helical_dee_source import RMIN,RMAX


def stencil(n,h,degree,periodic=False):
    count=5 if degree==1 or periodic else 6
    if n<count:raise ValueError('Grid too small for fixed fourth-order stencil')
    matrix=sp.lil_matrix((n,n))
    for i in range(n):
        start=min(max(i-count//2,0),n-count)
        js=np.arange(i-2,i+3) if periodic else np.arange(start,start+count)
        offsets=(js-i).astype(float);rhs=np.zeros(count);rhs[degree]=1 if degree==1 else 2
        weights=np.linalg.solve(np.array([offsets**k for k in range(count)]),rhs)/h**degree
        for j,w in zip(js,weights):matrix[i,j%n]+=w
    return matrix.tocsr()


class FourthVectorGrid(baseline.TensorGrid):
    def __init__(self,n,stabilization=.25):
        # The fixed baseline scalar flux still includes its original .25 term.
        super().__init__(n,stabilization)
        nr,nt,np_=self.shape;ir,it,ip=[sp.eye(k,format='csr') for k in self.shape]
        kron=lambda a,b,c:sp.kron(sp.kron(a,b,format='csr'),c,format='csr')
        hr=(RMAX-RMIN)/n;ha=2*np.pi/n
        dr=kron(stencil(nr,hr,1),it,ip);dt=kron(ir,stencil(nt,ha,1,True),ip);dp=kron(ir,it,stencil(np_,ha,1,True))
        drr=kron(stencil(nr,hr,2),it,ip);dtt=kron(ir,stencil(nt,ha,2,True),ip);dpp=kron(ir,it,stencil(np_,ha,2,True))
        self.d=[dr,dt,dp];r,t,q=self.r,self.t,self.q;diag=sp.diags
        er,et,ep=self.basis
        self.grad=[(diag(er[i])@dr+diag(et[i]/r)@dt+diag(ep[i]/q)@dp).tocsr() for i in range(3)]
        self.components=[drr,
            diag(1/r**2)@dtt+diag(1/r)@dr,
            diag(1/q**2)@dpp+diag(np.cos(t)/q)@dr-diag(np.sin(t)/(r*q))@dt,
            diag(1/r)@(dr@dt)-diag(1/r**2)@dt,
            diag(1/q)@(dr@dp)-diag(np.cos(t)/q**2)@dp,
            diag(1/(r*q))@(dt@dp)+diag(np.sin(t)/q**2)@dp]
        def apply(x):
            full=np.zeros((self.size,3));full[self.idx]=x.reshape(3,-1).T
            hc=[op@full for op in self.components]
            lap=hc[0]+hc[1]+hc[2]
            out=[]
            for i in range(3):
                value=lap[:,i].copy()
                for j in range(3):
                    coefficients=[er[i]*er[j],et[i]*et[j],ep[i]*ep[j],er[i]*et[j]+et[i]*er[j],er[i]*ep[j]+ep[i]*er[j],et[i]*ep[j]+ep[i]*et[j]]
                    value+=sum(c*h[:,j] for c,h in zip(coefficients,hc))/3
                out.append(value[self.idx])
            return np.concatenate(out)
        m=3*len(self.idx);self.vector=LinearOperator((m,m),matvec=apply,dtype=float)

    def solve_vector(self,u):
        w=super().solve_vector(u)
        rhs=(8*np.pi*baseline.G*u[self.idx]**6*self.S[:,self.idx]).ravel()
        residual=float(np.linalg.norm(self.vector@w[:,self.idx].ravel()-rhs));tol=max(1e-13,1e-10*float(np.linalg.norm(rhs)))
        if not np.isfinite(residual) or residual>tol:raise RuntimeError(f'True vector residual {residual} exceeds {tol}')
        self.history[-1]['acceptance_tolerance']=tol
        return w


if __name__=='__main__':
    original_solve=FourthVectorGrid.solve
    def solve_with_snapshot(self):
        out=baseline.output_dir()
        (out/'executed_fourth_vector.py').write_bytes(Path(__file__).read_bytes())
        (out/'effective-settings.json').write_text(json.dumps(dict(
            vector_order=4,vector_stabilization=0,
            scalar_order=2,scalar_stabilization=.25,
            radial_D1_points=5,radial_D2_points=6,
            periodic_D1_D2_points=5,
            note='Inherited settings.json stabilization=.25 describes the scalar operator; the direct vector operator has no stabilization.'),indent=2)+'\n')
        (out/'correction.json').write_text(json.dumps(dict(vector='fourth-order direct Cartesian Hessian with analytic chart connections',scalar='unchanged baseline stabilized second-order conservative face flux',physical_source='unchanged',frozen_audit='unchanged',vector_stabilization=0),indent=2)+'\n')
        return original_solve(self)
    FourthVectorGrid.solve=solve_with_snapshot
    baseline.TensorGrid=FourthVectorGrid
    baseline.main()
