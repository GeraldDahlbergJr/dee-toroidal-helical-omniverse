"""Independent analytic reference for the fourth-vector candidate's stencils."""
import sys,json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
from run43_fourth_vector import FourthVectorGrid
from diagnose_toroidal_boundary_convergence import Audit
from investigate_toroidal_inner_operator import region_norms
from validate_toroidal_helical_dee_source import RMIN,RMAX


def full_operator(g,w):
    # Evaluate stored chart Hessians on supplied analytic boundary values.
    hc=[op@w.T for op in g.components];er,et,ep=g.basis
    result=np.empty_like(w)
    for i in range(3):
        result[i]=(hc[0]+hc[1]+hc[2])[:,i]
        for j in range(3):
            result[i]+=(er[i]*er[j]*hc[0][:,j]+et[i]*et[j]*hc[1][:,j]+ep[i]*ep[j]*hc[2][:,j]+(er[i]*et[j]+et[i]*er[j])*hc[3][:,j]+(er[i]*ep[j]+ep[i]*er[j])*hc[4][:,j]+(et[i]*ep[j]+ep[i]*et[j])*hc[5][:,j])/3
    return result


def check(levels=(8,16,32,64)):
    rows=[]
    for n in levels:
        g=FourthVectorGrid(n);a=Audit(g.shape);r=a.r
        v=(r-RMIN)*(RMAX-r);vp=RMIN+RMAX-2*r;e=np.exp(8*(r-RMIN))
        f=v*e;fp=e*(vp+8*v);fpp=e*(-2+16*vp+64*v);c=np.array([.3,-.4,.5])[:,None,None,None]
        w=np.broadcast_to(c*f,(3,)+g.shape).copy();erdot=np.sum(a.er*c,axis=0);etdot=np.sum(a.et*c,axis=0);epdot=np.sum(a.ep*c,axis=0)
        truth=c*(fpp+(1/r+np.cos(a.t)/a.q)*fp)+(a.er*fpp*erdot+a.et*fp/r*etdot+a.ep*fp*np.cos(a.t)/a.q*epdot)/3
        actual=full_operator(g,w.reshape(3,-1)).reshape((3,)+g.shape)
        np.testing.assert_allclose((g.vector@w.reshape(3,-1)[:,g.idx].ravel()).reshape(3,-1),actual.reshape(3,-1)[:,g.idx],atol=1e-10,rtol=1e-10)
        norms=region_norms(a,actual-truth,np.ones(g.shape))
        norms['outer_first_row']=float(np.sqrt(np.mean((actual[:,-2]-truth[:,-2])**2)))
        norms['endpoint_max']=float(np.max(abs((actual-truth)[:,[0,-1]])))
        # Cartesian quadratic vector has exact Delta_L=(7,0,8/3), independent
        # of chart formula and with genuine angular dependence.
        x=g.q*np.cos(g.p);y=g.q*np.sin(g.p);z=g.r*np.sin(g.t)
        polynomial=np.array([x*x+y*y+z*z,x*y,z*z])
        error=full_operator(g,polynomial)-np.array([7.,0.,8/3])[:,None]
        rows.append(dict(n=n,radial_errors=norms,cartesian_polynomial_rms=float(np.sqrt(np.mean(error[:,g.idx]**2)))))
    orders={region:[float(np.log2(x['radial_errors'][region]/y['radial_errors'][region])) for x,y in zip(rows,rows[1:])] for region in rows[0]['radial_errors']}
    return dict(scope='Manufactured analytic operator accuracy, separate from physical-source validation.',rows=rows,radial_orders=orders,cartesian_polynomial_orders=[float(np.log2(x['cartesian_polynomial_rms']/y['cartesian_polynomial_rms'])) for x,y in zip(rows,rows[1:])])


if __name__=='__main__':
    report=check();Path(sys.argv[1]).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
