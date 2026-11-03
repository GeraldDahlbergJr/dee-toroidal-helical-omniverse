"""Executable discrete identity and smooth-field consistency; no physical solves.

Manufactured fields below test truncation only. They are never substituted for
original-source validation fields or for the frozen acceptance checks.
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'reproduce'))
from run43_consistent_momentum import TensorGrid
from run43_endpoint_continuation import EndpointGrid
from investigate_toroidal_inner_operator import tensor_div, region_norms
from diagnose_toroidal_boundary_convergence import Audit, radial_matrix
from validate_toroidal_helical_dee_source import RMIN, RMAX
from dee_stress_energy_3p1 import direct_tensor_check


def radial_errors():
    from investigate_toroidal_inner_operator import endpoint_cubic_matrix
    rows=[]
    for n in (8,16,32,64,128):
        r=np.linspace(RMIN,RMAX,n+1);h=(RMAX-RMIN)/n
        f=(r-RMIN)*(RMAX-r)*np.exp(8*(r-RMIN))
        v=(r-RMIN)*(RMAX-r);vp=RMIN+RMAX-2*r
        truth=np.exp(8*(r-RMIN))*(-2+16*vp+64*v)
        row=dict(n=n)
        for name,d in [('original',radial_matrix(n+1,h,2)),('matched',endpoint_cubic_matrix(n+1,h,True))]:
            error=d@d@f-truth
            row[name]=dict(first_inner=abs(float(error[1])),first_outer=abs(float(error[-2])),endpoint_max=float(max(abs(error[0]),abs(error[-1]))))
        rows.append(row)
    return rows


def main():
    rows=[]
    for n in (8,16,32):
        g=TensorGrid(n,stabilization=0);c=EndpointGrid(n,stabilization=0)
        a=Audit(g.shape);rng=np.random.default_rng(43)
        w=rng.normal(size=(3,g.size));w[:,np.setdiff1d(np.arange(g.size),g.idx)]=0
        dw=np.array([[g.grad[i]@w[j] for j in range(3)] for i in range(3)])
        div=sum(dw[i,i] for i in range(3))
        expanded=np.array([g.lap@w[i]+g.grad[i]@div/3 for i in range(3)])
        tensor=np.array([sum(g.grad[j]@g.longitudinal(w)[i,j] for j in range(3)) for i in range(3)])
        defect=np.array([sum(g.grad[j]@(g.grad[j]@w[i]) for j in range(3))-g.lap@w[i]
                         +sum(g.grad[j]@(g.grad[i]@w[j])-g.grad[i]@(g.grad[j]@w[j]) for j in range(3)) for i in range(3)])
        difference=(tensor-expanded)[:,g.idx];identity_error=float(np.max(abs(difference-defect[:,g.idx])))
        scale=float(np.max(abs(difference)))
        assert identity_error < 1e-10*max(1,scale)
        np.testing.assert_array_equal(c.S,g.S)
        u=1+1e-3*rng.random(g.size)
        np.testing.assert_array_equal(c.energy(u),g.energy(u))
        source_error=0.
        for k in (0,g.size//2,g.size-1):
            er,et,ep=g.basis
            from validate_toroidal_helical_dee_source import RHO,PI,KTH,KPS,M
            D=np.zeros((3,3));D[1]=KTH*(et[:,k]/g.r[k]-M*ep[:,k]/g.q[k]);D[2]=KPS*(et[:,k]/g.r[k]-M*ep[:,k]/g.q[k])
            energy,S,*_=direct_tensor_check(RHO,PI,D,np.eye(3)*u[k]**4)
            source_error=max(source_error,abs(energy-g.energy(u)[k]),float(np.max(abs(S-g.S[:,k]))))
        # Smooth radial Cartesian vector with exact Dirichlet values.
        r=a.r;v=(r-RMIN)*(RMAX-r);vp=RMIN+RMAX-2*r;e=np.exp(8*(r-RMIN))
        f=v*e;fp=e*(vp+8*v);fpp=e*(-2+16*vp+64*v)
        vector=np.array([.3,-.4,.5])[:,None,None,None]
        smooth=np.broadcast_to(vector*f,(3,)+a.shape).copy()
        erdot=np.sum(a.er*vector,axis=0);etdot=np.sum(a.et*vector,axis=0);epdot=np.sum(a.ep*vector,axis=0)
        truth=vector*(fpp+(1/r+np.cos(a.t)/a.q)*fp)+(a.er*fpp*erdot+a.et*fp/r*etdot+a.ep*fp*np.cos(a.t)/a.q*epdot)/3
        errors={}
        from investigate_toroidal_inner_operator import expanded as audit_expanded,endpoint_cubic_matrix
        errors['expanded']=region_norms(a,audit_expanded(a,smooth)-truth,np.ones(a.shape))
        errors['tensor']=region_norms(a,tensor_div(a,smooth)-truth,np.ones(a.shape))
        a.drad[2]=endpoint_cubic_matrix(n+1,a.dr,True)
        errors['matched_tensor']=region_norms(a,tensor_div(a,smooth)-truth,np.ones(a.shape))
        rows.append(dict(n=n,identity_max_error=identity_error,identity_relative_error=identity_error/max(1,scale),direct_4d_source_max_error=source_error,manufactured_errors=errors))
    report=dict(scope='Smooth manufactured operator tests only; original-source solves are recorded separately.',identity='T_i - E_i = (sum_j D_j D_j - L)W_i + sum_j [D_j,D_i]W_j',rows=rows,radial_composition=radial_errors())
    path=Path(sys.argv[1]);path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
