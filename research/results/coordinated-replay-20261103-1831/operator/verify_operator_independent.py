"""Independent algebra/continuum consistency checks; never physical acceptance solves."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source/reproduce'))
from diagnose_toroidal_boundary_convergence import Audit
FROZEN=ROOT/'source/reproduce/diagnose_toroidal_boundary_convergence.py'
b=FROZEN.read_bytes()
assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()=='5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
rows=[]
for n in [8,16,32,64]:
 a=Audit((n+1,n,n));r=a.r;h=.2/n;c=np.array([.3,-.4,.5])[:,None,None,None]
 # Analytic continuum W=c f(r), f=.3-r-.0125/r, zero at both boundaries.
 f=.3-r-.0125/r;fp=-1+.0125/r**2;fpp=-.025/r**3
 H=a.er[:,None]*a.er[None,:]*fpp+a.et[:,None]*a.et[None,:]*fp/r+a.ep[:,None]*a.ep[None,:]*fp*np.cos(a.t)/a.q
 forcing=c*np.trace(H,axis1=0,axis2=1)+np.einsum('ij...,j...->i...',H,c)/3
 w=np.broadcast_to(c*f,(3,)+a.shape).copy();t=a.tensor(w,2)
 T=np.array([sum(a.grad(t[i,j],2)[j] for j in range(3)) for i in range(3)])
 residual=T-forcing
 rv=a.r[:,0,0];band=np.broadcast_to(((rv>=.05-1e-12)&(rv<=.1+1e-12)&(rv>.05+1e-12))[:,None,None],a.shape)
 rms=a.norms(residual,np.ones(a.shape),band)['volume_weighted_rms']
 d=a.drad[2];d2=np.zeros_like(d)
 for i in range(1,n):d2[i,i-1:i+2]=[1,-2,1]
 d2/=h*h
 exact=np.zeros(n+1);exact[:4]=np.array([-1,3,-3,1])/(4*h*h)
 error=float(np.max(abs((d@d-d2)[1]-exact)))
 assert error<1e-8
 # Tensor-divergence identity with independently expanded commutator.
 dw=np.array([a.grad(v,2) for v in w]);div=np.trace(dw,axis1=0,axis2=1)
 composed=np.array([sum(a.grad(a.grad(v,2)[j],2)[j] for j in range(3)) for v in w])
 graddiv=a.grad(div,2)
 E=composed+graddiv/3
 comm=np.array([sum(a.grad(a.grad(w[j],2)[i],2)[j]-a.grad(a.grad(w[j],2)[j],2)[i] for j in range(3)) for i in range(3)])
 ident=float(np.max(abs(T-E-comm)))
 assert ident<1e-8
 rows.append(dict(n=n,first_row_coefficient_identity_error=error,tensor_commutator_identity_error=ident,exact_continuum_field_inner_frozen_M2_truncation_rms=rms))
 del H,T,w,t,residual,dw,composed,comm
orders=[float(np.log2(rows[i]['exact_continuum_field_inner_frozen_M2_truncation_rms']/rows[i+1]['exact_continuum_field_inner_frozen_M2_truncation_rms'])) for i in range(3)]
out=dict(scope='Manufactured analytic consistency only; not original-source physical solves or acceptance.',frozen_audit_git_blob='5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22',rows=rows,inner_orders=orders,threshold_unchanged=1.5,passed=sum(p>=1.5 for p in orders))
(ROOT/'independent-operator-verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
