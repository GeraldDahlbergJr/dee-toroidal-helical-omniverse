"""Array-only independent fourth Hessian applied to fresh physical solved fields."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source/reproduce'))
from diagnose_toroidal_boundary_convergence import Audit
source=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=False)
(out/'executed_verifier.py').write_bytes(Path(__file__).read_bytes())
def radial(n,h,degree):
 count=5 if degree==1 else 6;D=np.zeros((n,n))
 for k in range(n):
  start=min(max(k-count//2,0),n-count);js=np.arange(start,start+count);x=(js-k).astype(float)
  rhs=np.zeros(count);rhs[degree]=1 if degree==1 else 2
  D[k,js]=np.linalg.solve(np.array([x**j for j in range(count)]),rhs)/h**degree
 return D
rows=[]
for n in (8,16,32,64):
 p=source/f'fields_{n}.npz';raw=p.read_bytes()
 with np.load(p) as z:u=z['psi'];w=z['W_cartesian']
 a=Audit(u.shape);r,t,q=a.r,a.t,a.q;h=.2/n;da=2*np.pi/n
 dr=radial(n+1,h,1);drr=radial(n+1,h,2)
 R=lambda v:np.einsum('ij,jkl->ikl',dr,v)
 RR=lambda v:np.einsum('ij,jkl->ikl',drr,v)
 def D(v,axis):return (-np.roll(v,-2,axis)+8*np.roll(v,-1,axis)-8*np.roll(v,1,axis)+np.roll(v,2,axis))/(12*da)
 def DD(v,axis):return (-np.roll(v,-2,axis)+16*np.roll(v,-1,axis)-30*v+16*np.roll(v,1,axis)-np.roll(v,2,axis))/(12*da*da)
 er,et,ep=a.er,a.et,a.ep
 H=np.zeros_like(w)
 for j in range(3):
  v=w[j];vr=R(v);vt=D(v,1);vp=D(v,2)
  parts=[RR(v),DD(v,1)/r**2+vr/r,DD(v,2)/q**2+np.cos(t)*vr/q-np.sin(t)*vt/(r*q),R(vt)/r-vt/r**2,R(vp)/q-np.cos(t)*vp/q**2,D(vp,1)/(r*q)+np.sin(t)*vp/q**2]
  H[j]+=parts[0]+parts[1]+parts[2]
  for i in range(3):
   co=[er[i]*er[j],et[i]*et[j],ep[i]*ep[j],er[i]*et[j]+et[i]*er[j],er[i]*ep[j]+ep[i]*er[j],et[i]*ep[j]+ep[i]*et[j]]
   H[i]+=sum(c*part for c,part in zip(co,parts))/3
 tensor=a.tensor(w,2)
 T=np.array([sum(a.grad(tensor[i,j],2)[j] for j in range(3)) for i in range(3)])
 b=8*np.pi*.001*u**6*a.S
 algebraic=u**-10*(H-b);mismatch=u**-10*(T-H);frozen=a.momentum(u,w,2)
 closure=float(np.max(abs(frozen[:,1:-1]-algebraic[:,1:-1]-mismatch[:,1:-1])))
 assert closure<1e-12
 def rms(v):
  weight=np.broadcast_to(r*q,u.shape)[1]*u[1]**6
  return float(np.sqrt(np.sum(np.sum(v[:,1]**2,axis=0)*weight)/np.sum(weight)))
 row=dict(n=n,field_sha256=hashlib.sha256(raw).hexdigest(),identity_max_abs=closure,first_row_algebraic_physical_rms=rms(algebraic),first_row_operator_mismatch_rms=rms(mismatch),first_row_frozen_M2_rms=rms(frozen),independent_fourth_unweighted_algebraic_l2=float(np.linalg.norm((H-b)[:,1:-1])),first_row_radius=float(r[1,0,0]))
 rows.append(row);np.savez_compressed(out/f'row_{n}.npz',algebraic_physical=algebraic[:,1],operator_mismatch=mismatch[:,1],frozen=frozen[:,1]);print(json.dumps(row),flush=True)
report=dict(scope='Original-source fourth-vector saved fields at the explicit source path, independently reconstructed array fourth-Hessian operator; not new solutions or replacement acceptance.',source=str(source),rows=rows,identity='M_frozen=psi^-10(HW-b)+psi^-10(TW-HW)')
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
