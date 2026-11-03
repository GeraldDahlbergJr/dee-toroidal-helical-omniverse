"""Analytic operator checks only; manufactured fields are never gate solutions."""
import sys,json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
from run43_radial_consistency import RadialGrid
from run43_consistent_momentum import TensorGrid
from diagnose_toroidal_boundary_convergence import Audit
from investigate_toroidal_inner_operator import region_norms
from validate_toroidal_helical_dee_source import RMIN,RMAX


def check(levels=(8,16,32,64)):
    rows=[]
    for n in levels:
        g=RadialGrid(n,stabilization=0)
        a=Audit(g.shape);r=a.r;v=(r-RMIN)*(RMAX-r);vp=RMIN+RMAX-2*r;e=np.exp(8*(r-RMIN))
        f=v*e;fp=e*(vp+8*v);fpp=e*(-2+16*vp+64*v)
        c=np.array([.3,-.4,.5])[:,None,None,None]
        w=np.broadcast_to(c*f,(3,)+a.shape).copy()
        erdot=np.sum(a.er*c,axis=0);etdot=np.sum(a.et*c,axis=0);epdot=np.sum(a.ep*c,axis=0)
        truth=c*(fpp+(1/r+np.cos(a.t)/a.q)*fp)+(a.er*fpp*erdot+a.et*fp/r*etdot+a.ep*fp*np.cos(a.t)/a.q*epdot)/3
        full=w.reshape(3,-1);tensor=g.longitudinal(full)
        baseline=np.array([sum(g.grad[j]@tensor[i,j] for j in range(3)) for i in range(3)])
        defect=np.array([g.radial_defect@z for z in full]);er=g.basis[0]
        corrected=baseline+defect+er*np.sum(er*defect,axis=0)/3
        applied=(g.vector@full[:,g.idx].ravel()).reshape(3,-1)
        np.testing.assert_allclose(applied,corrected[:,g.idx],atol=1e-10,rtol=1e-10)
        norms={name:region_norms(a,val.reshape((3,)+g.shape)-truth,np.ones(g.shape)) for name,val in [('baseline',baseline),('corrected',corrected)]}
        norms['baseline']['outer_first_row']=float(np.sqrt(np.mean((baseline.reshape((3,)+g.shape)[:,-2]-truth[:,-2])**2)))
        norms['corrected']['outer_first_row']=float(np.sqrt(np.mean((corrected.reshape((3,)+g.shape)[:,-2]-truth[:,-2])**2)))
        rows.append(dict(n=n,errors=norms))
    orders={variant:{region:[float(np.log2(x['errors'][variant][region]/y['errors'][variant][region])) for x,y in zip(rows,rows[1:])] for region in rows[0]['errors'][variant]} for variant in ('baseline','corrected')}
    return dict(scope='Manufactured operator accuracy only, no original-source field replacement.',correction='fixed radial principal symbol times (Drr-Dr Dr)',rows=rows,orders=orders)


if __name__=='__main__':
    report=check();Path(sys.argv[1]).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
