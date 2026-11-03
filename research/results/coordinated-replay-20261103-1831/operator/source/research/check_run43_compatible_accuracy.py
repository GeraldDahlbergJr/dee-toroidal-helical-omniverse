"""Analytic boundary accuracy and source-preservation checks for the new closure."""
import json, os, sys, math
from pathlib import Path
import numpy as np
from run43_compatible_closure import CompatibleGrid, fourth_difference
from run43_consistent_momentum import TensorGrid
from check_run43_fourth_accuracy import full_operator
from diagnose_toroidal_boundary_convergence import Audit


def check(grid_type=CompatibleGrid):
    rows=[]
    for n in (8,16,32):
        g=grid_type(n); a=Audit(g.shape);r=a.r
        v=(r-.05)*(.25-r);vp=.30-2*r;e=np.exp(8*(r-.05))
        f=v*e;fp=e*(vp+8*v);fpp=e*(-2+16*vp+64*v)
        c=np.array([.3,-.4,.5])[:,None,None,None]
        w=np.broadcast_to(c*f,(3,)+g.shape).copy()
        truth=c*(fpp+(1/r+np.cos(a.t)/a.q)*fp)+(a.er*fpp*np.sum(a.er*c,axis=0)+a.et*fp/r*np.sum(a.et*c,axis=0)+a.ep*fp*np.cos(a.t)/a.q*np.sum(a.ep*c,axis=0))/3
        actual=(g.vector@w.reshape(3,-1)[:,g.idx].ravel()).reshape(3,n-1,n,n)
        error=actual-truth[:,1:-1]
        x=g.q*np.cos(g.p);y=g.q*np.sin(g.p);z=g.r*np.sin(g.t)
        polynomial=np.array([x*x+y*y+z*z,x*y,z*z])
        angular=full_operator(g,polynomial)
        if hasattr(g,'ko'):
            for axis,scale in enumerate((np.ones(g.size),g.r**-2,g.q**-2)):
                v=(g.ko[axis]@polynomial.T).T;eaxis=g.basis[axis]
                h=.2/n if axis==0 else 2*np.pi/n
                angular-=h*h/4*scale*(v+eaxis*np.sum(eaxis*v,axis=0)/3)
        angular-=np.array([7.,0.,8/3])[:,None]
        rows.append(dict(n=n,inner=float(np.sqrt(np.mean(error[:,0]**2))),outer=float(np.sqrt(np.mean(error[:,-1]**2))),rms=float(np.sqrt(np.mean(error**2))),angular_mixed=float(np.sqrt(np.mean(angular[:,g.idx]**2)))))
        if n==8:
            b=TensorGrid(n)
            np.testing.assert_array_equal(g.S,b.S)
            np.testing.assert_array_equal(g.energy(np.linspace(1,1.01,g.size)),b.energy(np.linspace(1,1.01,g.size)))
            assert (g.lap-b.lap).nnz==0
        d4=fourth_difference(n+1,1.)
        x=np.arange(n+1,dtype=float)
        for power in range(6):
            exact=np.zeros_like(x) if power<4 else math.factorial(power)/math.factorial(power-4)*x**(power-4)
            np.testing.assert_allclose(d4@x**power,exact,atol=3e-6,rtol=1e-7)
    orders={k:[float(np.log2(a[k]/b[k])) for a,b in zip(rows,rows[1:])] for k in ('inner','outer','rms','angular_mixed')}
    return dict(rows=rows,orders=orders,source_and_scalar_preserved=True)


if __name__=='__main__':
    result=check();Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
