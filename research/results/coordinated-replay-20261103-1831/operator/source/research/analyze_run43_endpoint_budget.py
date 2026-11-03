"""Saved fourth-vector fields: exact endpoint and supplementary error budgets.

No solves, audit edits, fitted corrections or changes to acceptance criteria.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,hashlib,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reproduce'),str(ROOT/'research')]
from diagnose_toroidal_boundary_convergence import Audit,radial_matrix
from investigate_toroidal_inner_operator import endpoint_cubic_matrix
from validate_original_physics import Independent


def derivative(n,h,degree,count):
    m=np.zeros((n,n))
    for i in range(n):
        js=np.arange(min(max(i-count//2,0),n-count),min(max(i-count//2,0),n-count)+count)
        x=(js-i).astype(float);b=np.zeros(count);b[degree]=math.factorial(degree)
        m[i,js]=np.linalg.solve(np.array([x**k for k in range(count)]),b)/h**degree
    return m


def momentum(a,u,w):
    dw=np.array([a.grad(v) for v in w]);div=np.trace(dw,axis1=0,axis2=1)
    t=dw+dw.swapaxes(0,1)-2/3*np.eye(3)[:,:,None,None,None]*div
    return u**-10*np.array([sum(a.grad(t[i,j])[j] for j in range(3)) for i in range(3)])-8*np.pi*.001*u**-4*a.S

class Hybrid:
    def __init__(self,a,ref):self.a=a;self.ref=ref;self.S=a.S
    def grad(self,v):
        a=self.a
        return a.er*self.ref.dr(v)+a.et*a.angular_derivative(v,1,a.dt,2)/a.r+a.ep*a.angular_derivative(v,2,a.dp,2)/a.q


def main():
    rows=[]
    for n in (8,16,32,64):
        p=ROOT/f'research/results/radial-retry-20261103-fourth-vector/fields_{n}.npz'
        with np.load(p) as z:u,w=z['psi'],z['W_cartesian']
        a=Audit(u.shape);ref=Independent(n);ref.S=a.S
        f=a.momentum(u,w,2)
        b=Audit(u.shape);b.drad[2]=endpoint_cubic_matrix(n+1,a.dr,True)
        matched=b.momentum(u,w,2);e=f-matched
        hybrid=momentum(Hybrid(a,ref),u,w);q=momentum(ref,u,w)
        terms={'endpoint':e,'remaining_radial':matched-hybrid,'angular':hybrid-q,'high_order_residual':q}
        weight=np.broadcast_to(a.r*a.q,u.shape)*u**6
        mask=np.broadcast_to((a.r>.05+1e-12)&(a.r<=.10+1e-12),u.shape)
        def dot(v,z):return float(np.sum((np.sum(v*z,axis=0)*weight)[mask])/np.sum(weight[mask]))
        def norm(v):return math.sqrt(max(0,dot(v,v)))
        d3=np.diff(w[:,:4],n=3,axis=1)[:,0]/a.dr**3
        exact=np.zeros_like(w);er=a.er[:,1]
        exact[:,1]=u[1]**-10*a.dr/4*(d3+er*np.sum(er*d3,axis=0)/3)
        for index in (1,-2):
            # Formula here asserts only the inner row; outer row excluded from inner band.
            if index==1:np.testing.assert_allclose(e[:,index],exact[:,index],atol=3e-14,rtol=1e-8)
        assert norm(e-exact)<1e-12
        predictions={}
        for count in (5,7,8):
            d=derivative(n+1,a.dr,3,count)
            d3hi=np.einsum('j,cjkl->ckl',d[1],w)
            v=np.zeros_like(w);v[:,1]=u[1]**-10*a.dr/4*(d3hi+er*np.sum(er*d3hi,axis=0)/3)
            predictions[str(count)]=dict(rms=norm(v),relative_error_vs_exact_endpoint=norm(v-e)/norm(e),cosine=dot(v,e)/norm(v)/norm(e))
        components={k:norm(v) for k,v in terms.items()}
        gram={k:{l:dot(v,z) for l,z in terms.items()} for k,v in terms.items()}
        total=sum(terms.values());np.testing.assert_allclose(total,f,atol=1e-15,rtol=1e-10)
        remainder=f-e;cross=2*dot(e,remainder);f2=dot(f,f)
        row1=np.zeros_like(w);row1[:,1]=f[:,1]
        row=dict(n=n,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),inner_samples=n//4,
                 frozen_rms=norm(f),matched_rms=norm(matched),hybrid_rms=norm(hybrid),high_order_rms=norm(q),
                 components=components,gram=gram,endpoint_cross_with_remainder=cross,
                 endpoint_remainder_cosine=dot(e,remainder)/norm(e)/norm(remainder),
                 cancellation_fraction=-cross/(dot(e,e)+dot(remainder,remainder)),
                 frozen_squared_identity=dict(endpoint=dot(e,e),remainder=dot(remainder,remainder),cross=cross,total=f2),
                 first_row_squared_fraction=dot(row1,row1)/f2,leading_predictions=predictions,
                 exact_endpoint_formula_relative_error=norm(e-exact)/norm(e))
        rows.append(row);print(n,row['frozen_rms'],row['components'],flush=True)
    pairs=[]
    for a,b in zip(rows,rows[1:]):
        order=lambda x,y:float(np.log2(x/y))
        ratio=lambda row:row['frozen_rms']/math.sqrt(row['frozen_squared_identity']['endpoint']+row['frozen_squared_identity']['remainder'])
        pairs.append(dict(pair=[a['n'],b['n']],frozen_order=order(a['frozen_rms'],b['frozen_rms']),
                          endpoint_order=order(a['components']['endpoint'],b['components']['endpoint']),
                          remainder_order=order(a['matched_rms'],b['matched_rms']),
                          incoherent_order=order(math.hypot(a['components']['endpoint'],a['matched_rms']),math.hypot(b['components']['endpoint'],b['matched_rms'])),
                          cancellation_order_shift=order(ratio(a),ratio(b)),
                          component_orders={k:order(a['components'][k],b['components'][k]) for k in a['components']}))
    out=dict(scope='Saved fourth-order vector candidate only; supplementary decomposition, unchanged frozen gate. Sequential split: endpoint, radial, angular, residual. Mixed effects depend on split order; exact endpoint identity does not.',rows=rows,pairs=pairs)
    target=ROOT/'research/results/endpoint-error-budget-20261103/budget.json'
    target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(pairs,indent=2))
if __name__=='__main__':main()
