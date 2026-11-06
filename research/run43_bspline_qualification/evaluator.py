"""Independent strong evaluator: piecewise-polynomial Horner + Cartesian chain rule.
Does not call Space.grad/scalar/vector or assembly quadrature."""
import numpy as np
from scipy.interpolate import BSpline,PPoly
from scipy.sparse import csr_matrix
from numpy.polynomial.legendre import leggauss
from experiment import A,B,exact_at,tensor

def polynomial_basis(n,x,axis,order):
    # Convert each scalar spline to local powers, then differentiate powers.
    if axis==0:
        knots=np.r_[[A]*5,np.linspace(A,B,n+1)[1:-1],[B]*5]; count=n+4
        mapping=list(range(1,count-1)); xx=np.asarray(x)
    else:
        knots=np.arange(-4,n+5)*(2*np.pi/n);count=n+4
        mapping=[(i-4)%n for i in range(count)];xx=np.asarray(x)%(2*np.pi)
    out=np.zeros((len(xx),n+2 if axis==0 else n))
    for i in range(count):
        if axis==0 and i not in mapping:continue
        coef=np.zeros(count);coef[i]=1
        pp=PPoly.from_spline((knots,coef,4)); ids=np.searchsorted(pp.x,xx,side='right')-1;ids=np.clip(ids,0,len(pp.x)-2)
        v=xx-pp.x[ids]; cc=pp.c.copy()
        for _ in range(order):cc=cc[:-1]*np.arange(len(cc)-1,0,-1)[:,None]
        value=np.zeros_like(xx)
        for row in cc:value=value*v+row[ids]
        out[:,i-1 if axis==0 else mapping[i]]+=value
    return csr_matrix(out)

DERS=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]

def geometry(r,t,p):
    q=1+r*np.cos(t);st=np.sin(t);ct=np.cos(t);sp=np.sin(p);cp=np.cos(p)
    er=np.array([ct*cp,ct*sp,st+0*p]); et=np.array([-st*cp,-st*sp,ct+0*p]);ep=np.array([-sp+0*t,cp+0*t,np.zeros_like(t+p)])
    shape=np.broadcast_shapes(r.shape,t.shape,p.shape)
    E=np.array([er,et/r,ep/q]) # coordinate a, Cartesian i
    dE=np.zeros((3,3,3)+shape) # derivative a, coordinate b, Cartesian i
    dE[1,0]=et;dE[2,0]=np.array([-ct*sp,ct*cp,0*t+0*p])
    dE[0,1]=-et/r**2;dE[1,1]=-er/r;dE[2,1]=np.array([st*sp,-st*cp,0*t+0*p])/r
    dE[0,2]=-ep*ct/q**2;dE[1,2]=ep*r*st/q**2;dE[2,2]=np.array([-cp+0*t,-sp+0*t,0*t+0*p])/q
    return E,dE,q

def fields(c,ms,r,t,p):
    ds=np.array([[tensor([ms[a][orders[a]] for a in range(3)],z)[0] for orders in DERS] for z in c])
    shape=ds.shape[2:];g=np.zeros((4,3)+shape);h=np.zeros((4,3,3)+shape)
    E,dE,q=geometry(r,t,p)
    first=ds[:,1:4];second=np.zeros((4,3,3)+shape)
    for a in range(3):second[:,a,a]=ds[:,4+a]
    for idx,(a,b) in enumerate([(0,1),(0,2),(1,2)]):second[:,a,b]=second[:,b,a]=ds[:,7+idx]
    for i in range(3):
        g[:,i]=sum(E[a,i]*first[:,a] for a in range(3))
        for j in range(3):
            h[:,i,j]=sum(E[a,j]*E[b,i]*second[:,a,b]+E[a,j]*dE[a,b,i]*first[:,b] for a in range(3) for b in range(3))
    u=1+ds[0,0];w=ds[1:,0]
    strain=g[1:]+g[1:].swapaxes(0,1);div=sum(g[i+1,i] for i in range(3))
    for i in range(3):strain[i,i]-=2*div/3
    lap=np.trace(h,axis1=1,axis2=2)
    vec=lap[1:]+np.array([sum(h[j+1,j,i] for j in range(3))/3 for i in range(3)])
    return u,w,strain,lap[0],vec,q

def evaluate(c,n,k):
    z,v=leggauss(k); rnodes=(np.linspace(A,B,n+1)[:-1,None]+(z+1)*(B-A)/(2*n)).ravel();rw=np.tile(v*(B-A)/(2*n),n)
    angles=(np.arange(n)[:,None]*2*np.pi/n+(z+1)*np.pi/n).ravel();aw=np.tile(v*np.pi/n,n)
    t=angles[:,None];p=angles[None,:]
    angular=[ [polynomial_basis(n,angles,a,o) for o in range(3)] for a in (1,2)]
    names=['inner','bulk','outer','first_inner','first_outer'];sums={key:np.zeros(4) for key in names};weights={key:0. for key in names};err=np.zeros(4);errw=0;maxerr=np.zeros(4)
    checks=[]
    for ir,rr in enumerate(np.r_[rnodes,A+(B-A)/n,B-(B-A)/n]):
        r=np.array([[rr]]);rad=[polynomial_basis(n,[rr],0,o) for o in range(3)]
        u,w,strain,lap,vec,q=fields(c,[rad]+angular,r,t,p)
        vals,gr,op=exact_at(r,t,p);du=gr[1:];div=sum(du[i,i] for i in range(3));es=du+du.swapaxes(0,1)
        for i in range(3):es[i,i]-=2*div/3
        bb=.068614375*(rr**-2+9/q**2)
        src=lambda uu:2*np.pi*.001*(.0127638*uu**5+bb*uu)
        E,_,_=geometry(r,t,p);S=-.030316*(E[1]-3*E[2])
        fh=-op[0]-np.sum(es**2,axis=(0,1))*vals[0]**-7/8-src(vals[0]);fm=-op[1:]+8*np.pi*.001*vals[0]**6*S
        H=-8*u**-5*(lap+np.sum(strain**2,axis=(0,1))*u**-7/8+src(u)+fh)
        M=u**-10*(vec-8*np.pi*.001*u**6*S+fm)
        res=np.concatenate([H[None],M]);ww=rr*q*u**6*aw[:,None]*aw[None,:]
        if ir<len(rnodes):
            keys=['inner' if rr<.10 else 'outer' if rr>.20 else 'bulk'];ww=ww*rw[ir]
            ee=np.concatenate([(u-vals[0])[None],w-vals[1:]])
            err+=np.sum(ee**2*ww,axis=(1,2));errw+=np.sum(ww);maxerr=np.maximum(maxerr,np.max(abs(ee),axis=(1,2)))
        else:keys=[names[ir-len(rnodes)+3]]
        for key in keys:sums[key]+=np.sum(res**2*ww,axis=(1,2));weights[key]+=np.sum(ww)
    return {'residual_rms':{key:np.sqrt(sums[key]/weights[key]).tolist() for key in names},'field_rms':np.sqrt(err/errw).tolist(),'field_max':maxerr.tolist()}

def traces(c,n):
    # Supremum certificate from endpoint/periodic coefficient construction;
    # evaluate a non-grid set too, including both seam copies.
    angles=np.linspace(0,2*np.pi,2*n+3)
    maxrad=np.zeros(4);maxseam=np.zeros(4);disagreement=0.
    from experiment import basis
    for axis in range(3):
        x=np.array([A,B]) if axis==0 else np.array([0,2*np.pi])
        mats=[polynomial_basis(n,x if a==axis else (np.linspace(A,B,7) if a==0 else angles),a,0) for a in range(3)]
        for i,z in enumerate(c):
            val=tensor(mats,z)
            if axis==0:maxrad[i]=max(maxrad[i],np.max(abs(val)))
            else:maxseam[i]=max(maxseam[i],np.max(abs(np.take(val,0,axis=axis)-np.take(val,1,axis=axis))))
    for a in range(3):
        x=np.linspace(A,B,2*n+3) if a==0 else angles
        for o in range(3):
            disagreement=max(disagreement,np.max(abs((polynomial_basis(n,x,a,o)-basis(n,x,a,o)).toarray())))
    return {'radial_max':maxrad.tolist(),'seam_max':maxseam.tolist(),'basis_derivative_disagreement_max':float(disagreement),'scales':[1.,.0025,.00175,.00225],'tolerances':(100*np.finfo(float).eps*np.array([1.,.0025,.00175,.00225])).tolist()}
