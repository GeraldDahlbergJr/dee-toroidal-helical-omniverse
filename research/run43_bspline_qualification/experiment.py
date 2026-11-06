"""Single Run43 manufactured qualification; no baseline solver imports."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import json,time,resource,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator,cg
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
A=.05; B=.25; P=4; G=.001
START=time.monotonic()

def basis(n,x,axis,der=0):
    lo,hi=(A,B) if axis==0 else (0,2*np.pi)
    if axis==0:
        knots=np.r_[[lo]*5,np.linspace(lo,hi,n+1)[1:-1],[hi]*5]
        c=np.eye(n+4)
        return csr_matrix(BSpline(knots,c,4)(x,nu=der)[:,1:-1])
    h=(hi-lo)/n
    knots=np.arange(-4,n+5)*h
    c=np.zeros((n+4,n))
    for i in range(n+4): c[i,(i-4)%n]=1
    return csr_matrix(BSpline(knots,c,4)(np.asarray(x)%(2*np.pi),nu=der))

def rule(n,k,axis):
    l,h=(A,B) if axis==0 else (0,2*np.pi)
    z,w=leggauss(k); edges=np.linspace(l,h,n+1)
    return ((edges[:-1,None]+(z+1)*(h-l)/(2*n)).ravel(),np.tile(w*(h-l)/(2*n),n))

def apply_axis(m,x,a):
    y=np.moveaxis(x,a,0); shape=y.shape
    y=m@y.reshape(shape[0],-1)
    return np.moveaxis(np.asarray(y).reshape((m.shape[0],)+shape[1:]),0,a)

def tensor(ms,x,transpose=False):
    for a,m in enumerate(ms): x=apply_axis(m.T if transpose else m,x,a)
    return x

class Space:
    def __init__(self,n,k):
        self.n=n; self.k=k; self.cshape=(n+2,n,n); self.size=np.prod(self.cshape)
        rules=[rule(n,k,a) for a in range(3)]
        self.x=[z[0] for z in rules]
        r,t,p=np.meshgrid(*self.x,indexing='ij',sparse=True)
        self.r=r; self.t=t; self.p=p; self.q=1+r*np.cos(t)
        self.weight=rules[0][1][:,None,None]*rules[1][1][None,:,None]*rules[2][1][None,None,:]*r*self.q
        self.ms=[basis(n,self.x[a],a) for a in range(3)]
        self.ds=[basis(n,self.x[a],a,1) for a in range(3)]
        self.frame=[(np.cos(t)*np.cos(p),np.cos(t)*np.sin(p),np.sin(t)+0*p),(-np.sin(t)*np.cos(p)/r,-np.sin(t)*np.sin(p)/r,np.cos(t)/r+0*p),(-np.sin(p)/self.q,np.cos(p)/self.q,np.zeros_like(r+t+p))]
        self.prebuild(rules)
    def val(self,c): return tensor(self.ms,c)
    def derivatives(self,c):
        return [tensor([self.ds[a] if j==a else self.ms[j] for j in range(3)],c) for a in range(3)]
    def grad(self,c):
        d=self.derivatives(c)
        return [sum(self.frame[a][j]*d[a] for a in range(3)) for j in range(3)]
    def back(self,v): return tensor(self.ms,self.weight*v,True)
    def adj(self,v):
        out=np.zeros(self.cshape)
        for a in range(3):
            z=sum(self.frame[a][j]*v[j] for j in range(3))*self.weight
            out+=tensor([self.ds[a] if j==a else self.ms[j] for j in range(3)],z,True)
        return out
    def scalar(self,c): return self.adj(self.grad(c))
    def vector(self,c):
        # Direct continuum L: Cartesian symmetric trace-free strain.
        dw=np.array([self.grad(z) for z in c]); div=sum(dw[i,i] for i in range(3))
        strain=dw+dw.swapaxes(0,1)
        for i in range(3): strain[i,i]-=2*div/3
        return np.array([self.adj(strain[i]) for i in range(3)])
    def prebuild(self,rules):
        self.evals=[]; self.evec=[]
        for a in range(3):
            w=rules[a][1]; m=self.ms[a]; d=self.ds[a]
            mass=(m.T@ m.multiply(w[:,None])).toarray()
            stiff=(d.T@d.multiply(w[:,None])).toarray()
            ev,v=eigh(stiff,mass); self.evals.append(ev);self.evec.append(v)
        self.denom=self.evals[0][:,None,None]+self.evals[1][None,:,None]/.15**2+self.evals[2][None,None,:]
    def pre(self,c):
        z=tensor([csr_matrix(v.T) for v in self.evec],c)
        return tensor([csr_matrix(v) for v in self.evec],z/self.denom)/.15
    def sources(self,u):
        e0=.0127638; bb=.068614375*(self.r**-2+9*self.q**-2)
        S=np.array([-.030316*(self.frame[1][j]-3*self.frame[2][j])+np.zeros_like(self.weight) for j in range(3)])
        return 2*np.pi*G*(e0*u**5+bb*u),8*np.pi*G*u**6*S

def analytic_functions():
    # Symbolic continuum forcing, independent of Galerkin/basis differentiation.
    import sympy as s
    r,t,p=s.symbols('r t p',real=True); q=1+r*s.cos(t)
    f=s.sin(s.pi*(r-s.Rational(1,20))/s.Rational(1,5))*s.exp(s.Rational(1,5)*s.sin(t)+s.Rational(1,7)*s.cos(p))
    fields=[1+s.Rational(1,100)*f,s.Rational(1,1000)*f*(1+s.sin(t+p)),s.Rational(7,10000)*f*(1+s.cos(t-2*p)),s.Rational(9,10000)*f*(1+s.sin(2*t-p))]
    er=[s.cos(t)*s.cos(p),s.cos(t)*s.sin(p),s.sin(t)]
    et=[-s.sin(t)*s.cos(p)/r,-s.sin(t)*s.sin(p)/r,s.cos(t)/r]
    ep=[-s.sin(p)/q,s.cos(p)/q,0]
    def D(v,j): return er[j]*s.diff(v,r)+et[j]*s.diff(v,t)+ep[j]*s.diff(v,p)
    du=[[D(v,j) for j in range(3)] for v in fields]
    lap=lambda v:s.diff(v,r,2)+(1/r+s.cos(t)/q)*s.diff(v,r)+s.diff(v,t,2)/r**2-s.sin(t)*s.diff(v,t)/(r*q)+s.diff(v,p,2)/q**2
    div=sum(du[i+1][i] for i in range(3))
    op=[lap(fields[0])]+[lap(fields[i+1])+D(div,i)/3 for i in range(3)]
    deriv=[[s.diff(v,r,aa,t,bb,p,cc) for aa,bb,cc in [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]] for v in fields]
    funcs=[s.lambdify((r,t,p),expr,'numpy',cse=True) for expr in (fields,du,op)]
    return funcs

FUN=analytic_functions()

def exact_at(r,t,p):
    shape=np.broadcast_shapes(r.shape,t.shape,p.shape)
    vals=np.array([np.broadcast_to(z,shape) for z in FUN[0](r,t,p)])
    grad=np.array([[np.broadcast_to(z,shape) for z in row] for row in FUN[1](r,t,p)])
    op=np.array([np.broadcast_to(z,shape) for z in FUN[2](r,t,p)])
    return vals,grad,op

def forcing(space):
    vals,gr,op=exact_at(space.r,space.t,space.p); u=vals[0]
    dw=gr[1:]; div=sum(dw[i,i] for i in range(3)); strain=dw+dw.swapaxes(0,1)
    for i in range(3):strain[i,i]-=2*div/3
    sc,vc=space.sources(u)
    return -op[0]-np.sum(strain**2,axis=(0,1))*u**-7/8-sc, -op[1:]+vc

def linear(space,fn,rhs,vector=False,potential=None):
    shape=((3,)+space.cshape) if vector else space.cshape; N=int(np.prod(shape))
    op=LinearOperator((N,N),matvec=lambda x:fn(x.reshape(shape)).ravel())
    pre=LinearOperator((N,N),matvec=lambda x:(np.array([space.pre(z) for z in x.reshape(shape)]) if vector else space.pre(x.reshape(shape))).ravel())
    count=[0]
    def cb(x):count[0]+=1
    sol,info=cg(op,rhs.ravel(),M=pre,rtol=2e-12,atol=0,maxiter=500,callback=cb)
    rel=float(np.linalg.norm(op@sol-rhs.ravel())/np.linalg.norm(rhs))
    print(json.dumps({'linear': 'vector' if vector else 'scalar','n':space.n,'iterations':count[0],'info':info,'relative':rel,'elapsed':time.monotonic()-START}),flush=True)
    if info or rel>1e-10: raise RuntimeError('decisive linear algebraic failure')
    return sol.reshape(shape),rel

def solve(space):
    fh,fm=forcing(space); c=np.zeros(space.cshape); w=np.zeros((3,)+space.cshape); residuals=[]
    for it in range(12):
        u=1+space.val(c); sc,vc=space.sources(u)
        w,rel=linear(space,space.vector,np.array([space.back(fm[i]-vc[i]) for i in range(3)]),True);residuals.append(rel)
        dw=np.array([space.grad(z) for z in w]); div=sum(dw[i,i] for i in range(3)); strain=dw+dw.swapaxes(0,1)
        for i in range(3):strain[i,i]-=2*div/3
        sq=np.sum(strain**2,axis=(0,1))
        for j in range(8):
            u=1+space.val(c); sc,_=space.sources(u)
            rhs=space.back(fh+sq*u**-7/8+sc)-space.scalar(c)
            bb=.068614375*(space.r**-2+9*space.q**-2)
            pot=-7*sq*u**-8/8+2*np.pi*G*(5*.0127638*u**4+bb)
            step,rel=linear(space,lambda z:space.scalar(z)-space.back(pot*space.val(z)),rhs);residuals.append(rel); c+=step
            if np.max(abs(step))<2e-13: break
        change=np.max(abs(c-(old if it else np.zeros_like(c))))
        print(json.dumps({'coupled':it,'n':space.n,'change':float(change)}),flush=True)
        if it and change<2e-13:break
        old=c.copy()
    else:raise RuntimeError('coupled iteration failed')
    return np.concatenate([c[None],w]),residuals

if __name__=='__main__':
    print('implementation module; use run.py',flush=True)
