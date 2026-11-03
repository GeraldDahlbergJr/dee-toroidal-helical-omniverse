"""Factorial gradient-error expansion on saved fourth-vector fields only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,math,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'research'),str(ROOT/'reproduce')]
from analyze_run43_endpoint_budget import derivative
from diagnose_toroidal_boundary_convergence import Audit
from investigate_toroidal_inner_operator import endpoint_cubic_matrix
from validate_original_physics import Independent

def tensor(g,w):
    dw=np.array([g(v) for v in w]);div=np.trace(dw,axis1=0,axis2=1)
    return dw+dw.swapaxes(0,1)-2/3*np.eye(3)[:,:,None,None,None]*div

def divergence(g,t):return np.array([sum(g(t[i,j])[j] for j in range(3)) for i in range(3)])

def main():
    rows=[]
    for n in (8,16,32,64):
        path=ROOT/f'research/results/radial-retry-20261103-fourth-vector/fields_{n}.npz'
        with np.load(path) as z:u,w=z['psi'],z['W_cartesian']
        a=Audit(u.shape);ref=Independent(n);h=a.dr
        matched=endpoint_cubic_matrix(n+1,h,True)
        drdiff=matched-ref.d1
        gr=lambda v:a.er*np.einsum('ij,jkl->ikl',drdiff,v)
        ga=lambda v:a.et*(a.angular_derivative(v,1,a.dt,2)-ref.angle(v,1))/a.r+a.ep*(a.angular_derivative(v,2,a.dp,2)-ref.angle(v,2))/a.q
        g=ref.grad;t=tensor(g,w);tr=tensor(gr,w);ta=tensor(ga,w)
        linear=lambda delta,td:u**-10*(divergence(g,td)+divergence(delta,t))
        terms=dict(radial_linear=linear(gr,tr),radial_quadratic=u**-10*divergence(gr,tr),angular_linear=linear(ga,ta),angular_quadratic=u**-10*divergence(ga,ta),radial_angular_mixed=u**-10*(divergence(gr,ta)+divergence(ga,tr)))
        m=Audit(u.shape);m.drad[2]=matched
        q=u**-10*divergence(g,t)-8*np.pi*.001*u**-4*a.S
        target=m.momentum(u,w,2)-q
        weight=np.broadcast_to(a.r*a.q,u.shape)*u**6
        mask=np.broadcast_to((a.r>.05+1e-12)&(a.r<=.10+1e-12),u.shape)
        def dot(x,y):return float(np.sum((np.sum(x*y,axis=0)*weight)[mask])/np.sum(weight[mask]))
        norm=lambda x:math.sqrt(max(0,dot(x,x)))
        closure=norm(sum(terms.values())-target)
        assert closure<1e-12
        d3=derivative(n+1,h,3,7)
        grlead=lambda v:a.er*(h*h/6)*np.einsum('ij,jkl->ikl',d3,v)
        galead=lambda v:a.et*(a.dt*a.dt/6)*ref.angle(v,1,3)/a.r+a.ep*(a.dp*a.dp/6)*ref.angle(v,2,3)/a.q
        leading={}
        for name,delta in [('radial',grlead),('angular',galead)]:
            prediction=linear(delta,tensor(delta,w));actual=terms[name+'_linear']
            leading[name]=dict(predicted_rms=norm(prediction),actual_rms=norm(actual),relative_difference=norm(prediction-actual)/norm(actual),cosine=dot(prediction,actual)/norm(prediction)/norm(actual))
        # Taylor series of the exact endpoint delta3, evaluated at row 1.
        er=a.er[:,1];d3exact=np.diff(w[:,:4],n=3,axis=1)[:,0]/h**3
        project=lambda v:u[1]**-10*h/4*(v+er*np.sum(er*v,axis=0)/3)
        exact=np.zeros_like(w);exact[:,1]=project(d3exact)
        expansions={}
        for count in (7,8):
            cumulative=np.zeros_like(d3exact);steps=[]
            for k in range(3,count):
                coeff=float(np.dot(np.array([-1,3,-3,1]),np.array([-1.,0,1,2])**k)/math.factorial(k))
                dk=derivative(n+1,h,k,count)
                contribution=coeff*h**(k-3)*np.einsum('j,cjkl->ckl',dk[1],w)
                cumulative+=contribution;v=np.zeros_like(w);v[:,1]=project(cumulative)
                steps.append(dict(through_derivative=k,delta3_coefficient=coeff,relative_error=norm(v-exact)/norm(exact),rms=norm(v)))
            expansions[str(count)]=steps
        profiles={}
        for name,v in terms.items():
            square=np.sum(v*v,axis=0)*weight
            energy=np.sum(square[mask])
            profiles[name]=[dict(row=i,r=float(a.r[i,0,0]),squared_norm_fraction=float(np.sum(square[i])/energy)) for i in range(1,n//4+1)]
        row=dict(n=n,field_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),matched_minus_reference_rms=norm(target),components={k:norm(v) for k,v in terms.items()},gram={k:{l:dot(v,z) for l,z in terms.items()} for k,v in terms.items()},leading_predictions=leading,endpoint_taylor=expansions,identity_absolute_error=closure,radial_profiles=profiles)
        rows.append(row);print(n,json.dumps({k:row[k] for k in ('components','leading_predictions','endpoint_taylor')}),flush=True)
    report=dict(scope='Factorial radial/angular gradient perturbations relative to seven-point/Fourier reference. Supplementary only; no solver or acceptance changes. Taylor completion is interpolation identity, not continuum validation.',rows=rows)
    out=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'research/results/finite-grid-terms-20261103'
    out.mkdir(exist_ok=False)
    (out/'terms.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
