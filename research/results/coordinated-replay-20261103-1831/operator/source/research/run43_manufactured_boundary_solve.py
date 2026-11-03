"""New verification: absolute solution error versus frozen audit truncation.

Analytic momentum-only test problems, NEVER original-source acceptance fields.
No production source array or solver source is modified. The two forcing
functions are analytic Delta_L(c F), not fitted discrete forcing functions.
"""
import sys,json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import gmres,LinearOperator
from run43_consistent_momentum import TensorGrid
from run43_fourth_vector import FourthVectorGrid
from diagnose_toroidal_boundary_convergence import Audit


def analytic(g,kind):
    r=g.r;a=.05;b=.25;c=np.array([.3,-.4,.5])
    er,et,ep=g.basis;ct=np.cos(g.t);st=np.sin(g.t)
    f=a+b-r-a*b/r;fp=-1+a*b/r**2;fpp=-2*a*b/r**3
    # Independent Cartesian chain rule: Hess r = et et/r + cos(t) ep ep/q.
    hessr=np.einsum('in,jn->ijn',et,et)/r+ct/g.q*np.einsum('in,jn->ijn',ep,ep)
    if kind=='radial':
        F=f
        hf=fpp*np.einsum('in,jn->ijn',er,er)+fp*hessr
    else:
        # F=z*g(r), g=f/r. Differentiating in Cartesian coordinates gives
        # Hess F = g'(er ez^T+ez er^T)+z(g'' er er^T+g' Hess r).
        z=r*st;v=f/r;vp=-(a+b)/r**2+2*a*b/r**3;vpp=2*(a+b)/r**3-6*a*b/r**4
        ez=np.array([0.,0.,1.])
        F=z*v
        hf=vp*(np.einsum('in,j->ijn',er,ez)+np.einsum('i,jn->ijn',ez,er))+z*(vpp*np.einsum('in,jn->ijn',er,er)+vp*hessr)
    truth=c[:,None]*np.trace(hf,axis1=0,axis2=1)+np.einsum('ijn,j->in',hf,c)/3
    exact=c[:,None]*F
    # Floating subtraction in the closed form must not alter Dirichlet data.
    exact.reshape((3,)+g.shape)[:,[0,-1]]=0
    return exact,truth


def norms(a,v,band):
    rv=a.r[:,0,0];lo,hi=(.05,.10) if band=='inner' else (.20,.25)
    mask=np.broadcast_to(((rv>=lo-1e-12)&(rv<=hi+1e-12))[:,None,None],a.shape).copy();mask[[0,-1]]=False
    return a.norms(v.reshape((3,)+a.shape),np.ones(a.shape),mask)['volume_weighted_rms']


def main():
    out=Path(sys.argv[1])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    (out/'executed_manufactured_solve.py').write_bytes(Path(__file__).read_bytes())
    rows=[]
    for kind in ('radial','mixed'):
        for n in (8,16,32,64):
            for name,cls in [('tensor',TensorGrid),('fourth',FourthVectorGrid)]:
                start=time.monotonic();g=cls(n);a=Audit(g.shape)
                exact,truth=analytic(g,kind);rhs=truth[:,g.idx].ravel()
                ni=len(g.idx);pre=LinearOperator(g.vector.shape,matvec=lambda x:np.concatenate([g.ilu.solve(y) for y in x.reshape(3,ni)]),dtype=float)
                hist=[]
                x,info=gmres(g.vector,rhs,M=pre,rtol=1e-11,atol=1e-13,restart=60,maxiter=100,callback=lambda r:hist.append(float(r)),callback_type='pr_norm')
                residual=float(np.linalg.norm(g.vector@x-rhs));tol=max(1e-13,1e-11*np.linalg.norm(rhs))
                assert info==0 and np.isfinite(residual) and residual<=tol,(kind,n,name,info,residual,tol)
                solved=np.zeros_like(exact);solved[:,g.idx]=x.reshape(3,-1)
                err=solved-exact
                def tensor_div(w):
                    # Audit momentum subtracts its original matter source.
                    # Add it back ONLY to obtain the pure LINEAR test operator.
                    return (a.momentum(np.ones(g.shape),w.reshape((3,)+g.shape),2)+8*np.pi*.001*a.S).reshape(3,-1)
                # Remove any roundoff source cancellation by independently
                # forming div(LW) directly from the frozen tensor and gradients.
                def pure(w):
                    t=a.tensor(w.reshape((3,)+g.shape),2)
                    return np.array([sum(a.grad(t[i,j],2)[j] for j in range(3)) for i in range(3)]).reshape(3,-1)
                audit_consistency=pure(exact)-truth
                solution_effect=pure(err)
                measured=pure(solved)-truth
                closure=float(np.max(abs((measured-audit_consistency-solution_effect)[:,g.idx])))
                assert closure<1e-8
                bands={}
                for band in ('inner','outer'):
                    bands[band]={k:norms(a,v,band) for k,v in [('solution_error',err),('audit_on_exact',audit_consistency),('audit_on_solved',measured),('solution_effect',solution_effect)]}
                path=out/f'{kind}_{name}_{n}.npz'
                np.savez_compressed(path,exact=exact.reshape((3,)+g.shape),solved=solved.reshape((3,)+g.shape),analytic_forcing=truth.reshape((3,)+g.shape))
                row=dict(kind=kind,solver=name,n=n,bands=bands,max_solution_error=float(np.max(abs(err))),linear_residual=residual,tolerance=float(tol),gmres=hist,identity_max_abs=closure,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.monotonic()-start)
                rows.append(row);(out/f'{kind}_{name}_{n}.json').write_text(json.dumps(row,indent=2)+'\n')
                print(json.dumps({k:v for k,v in row.items() if k!='gmres'}),flush=True)
    orders=[]
    for kind in ('radial','mixed'):
        for name in ('tensor','fourth'):
            subset=[x for x in rows if x['kind']==kind and x['solver']==name]
            for band in ('inner','outer'):
                for metric in subset[0]['bands'][band]:
                    orders.append(dict(kind=kind,solver=name,band=band,metric=metric,orders=[float(np.log2(a['bands'][band][metric]/b['bands'][band][metric])) for a,b in zip(subset,subset[1:])]))
    (out/'report.json').write_text(json.dumps(dict(scope='Manufactured analytic momentum boundary-value problems only; not original-source acceptance or coupled production solves.',rows=rows,orders=orders),indent=2)+'\n')


if __name__=='__main__':main()
