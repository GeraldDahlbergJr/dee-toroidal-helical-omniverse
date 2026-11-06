"""One non-polynomial qualification; no solver imports or saved-field reads.
The direct-matrix ball calculation bounds the exact discrete operator; comparing
its center with the FFT evaluator includes observed FFT/library discrepancies.
"""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'run43_final_validation'))
from evaluate import Evaluator
OUT=Path(__file__).resolve().parent
mp.mp.dps=100
I=mp.iv
I.dps=80
EPS=np.finfo(float).eps

# Positive bound arithmetic uses at most 65-term contractions and fewer than
# 32 additional scalar operations per update. 1024 eps inflation dominates
# their gamma bounds and nextafter ensures outward final rounding.
def up(x): return np.nextafter(np.asarray(x)*(1+1024*EPS)+1e-290,np.inf)
class B:
    def __init__(self,c,e=0): self.c=np.asarray(c,dtype=float);self.e=np.broadcast_to(np.asarray(e,dtype=float),self.c.shape)
    def __add__(self,o):
        o=ball(o);c=self.c+o.c
        return B(c,up(self.e+o.e+EPS*abs(c)))
    __radd__=__add__
    def __neg__(self):return B(-self.c,self.e)
    def __sub__(self,o):return self+-ball(o)
    def __rsub__(self,o):return ball(o)+-self
    def __mul__(self,o):
        o=ball(o);c=self.c*o.c
        return B(c,up(abs(self.c)*o.e+abs(o.c)*self.e+self.e*o.e+EPS*abs(c)))
    __rmul__=__mul__
    def inv(self):
        lo=abs(self.c)-self.e
        assert np.all(lo>0)
        c=1/self.c
        return B(c,up(self.e/(abs(self.c)*lo)+EPS*abs(c)))
    def __truediv__(self,o):return self*ball(o).inv()
    def __rtruediv__(self,o):return ball(o)*self.inv()
    def __pow__(self,k):
        if k<0:return self.inv()**(-k)
        out=B(1)
        for _ in range(k):out=out*self
        return out

def ball(x):return x if isinstance(x,B) else B(x)
def exact(x):
    # Interval seeds, including transcendental functions, evaluated with
    # mpmath.iv at 80 decimal digits. Enclose endpoints in a binary64 ball.
    if not hasattr(x,'_mpi_'): x=I.mpf(x)
    lo=mp.make_mpf(x._mpi_[0]);hi=mp.make_mpf(x._mpi_[1])
    c=float((lo+hi)/2)
    return B(c,up(float(max(abs(mp.mpf(c)-lo),abs(mp.mpf(c)-hi)))))
def mp_array(values,shape):
    bs=[exact(v) for v in values]
    return B(np.array([v.c for v in bs]).reshape(shape),np.array([v.e for v in bs]).reshape(shape))
def linear(a,v,axis):
    a=ball(a);v=ball(v)
    c=np.moveaxis(np.tensordot(a.c,v.c,axes=(1,axis)),0,axis)
    mag=np.moveaxis(np.tensordot(abs(a.c),abs(v.c),axes=(1,axis)),0,axis)
    err=np.moveaxis(np.tensordot(abs(a.c),v.e,axes=(1,axis))+np.tensordot(a.e,abs(v.c)+v.e,axes=(1,axis)),0,axis)
    gamma=a.c.shape[1]*EPS/(1-a.c.shape[1]*EPS)
    return B(c,up(err+gamma*mag))

def setup(n):
    r=mp_array([I.mpf(1)/20+I.mpf(i)/(5*n) for i in range(n+1)],(n+1,1,1))
    trig={}
    for axis,shape in [(1,(1,n,1)),(2,(1,1,n))]:
        trig[axis]=[mp_array([f(2*I.pi*i/n) for i in range(n)],shape) for f in (I.sin,I.cos)]
    st,ct=trig[1];sp,cp=trig[2];q=1+r*ct
    er=[ct*cp,ct*sp,st];et=[-st*cp,-st*sp,ct];ep=[-sp,cp,B(0)]
    exp=mp_array([I.exp(12*(I.mpf(i)/n-1)) for i in range(n+1)],(n+1,1,1))
    a=1+exact(I.mpf(1)/5)*st+exact(I.mpf(3)/20)*cp+exact(I.mpf(1)/10)*(st*cp+ct*sp)
    at=exact(I.mpf(1)/5)*ct+exact(I.mpf(1)/10)*(ct*cp-st*sp)
    ap=-exact(I.mpf(3)/20)*sp+exact(I.mpf(1)/10)*(ct*cp-st*sp)
    att=-exact(I.mpf(1)/5)*st-exact(I.mpf(1)/10)*(st*cp+ct*sp)
    app=-exact(I.mpf(3)/20)*cp-exact(I.mpf(1)/10)*(st*cp+ct*sp)
    atp=-exact(I.mpf(1)/10)*(st*cp+ct*sp)
    f=exp*a;ft=exp*at;fp=exp*ap
    jet=[f,60*f,ft,fp,3600*f,60*ft,60*fp,exp*att,exp*atp,exp*app]
    # Reference derivative matrices use exact rational radial coefficients and
    # analytic real Fourier collocation matrices (even n, zero D1 Nyquist).
    rad=[]
    import sympy as sy
    for d in (1,2):
        mat=np.zeros((n+1,n+1));err=mat.copy()
        for i in range(n+1):
            start=min(max(i-4,0),n-8);ids=list(range(start,start+9))
            weights=sy.finite_diff_weights(d,[k-i for k in ids],0)[d][-1]
            for k,w in zip(ids,weights):
                rat=w*(5*n)**d;v=exact(I.mpf(int(sy.numer(rat)))/int(sy.denom(rat)))
                mat[i,k]=v.c;err[i,k]=v.e
        rad.append(B(mat,err))
    angular=[]
    for d in (1,2):
        mat=np.zeros((n,n));err=mat.copy()
        for j in range(n):
            for k in range(n):
                delta=j-k
                val=(I.mpf(0) if d==1 else -I.mpf(n*n+2)/12) if delta==0 else ((I.mpf(-1)**delta/2)*I.cos(I.pi*delta/n)/I.sin(I.pi*delta/n) if d==1 else -I.mpf(-1)**delta/(2*I.sin(I.pi*delta/n)**2))
                v=exact(val);mat[j,k]=v.c;err[j,k]=v.e
        angular.append(B(mat,err))
    return r,st,ct,q,er,et,ep,jet,rad,angular

def calculate(n):
    r,st,ct,q,er,et,ep,j,rad,ang=setup(n)
    def grad(v):return [er[i]*linear(rad[0],v,0)+et[i]*linear(ang[0],v,1)/r+ep[i]*linear(ang[0],v,2)/q for i in range(3)]
    def lap(v):return linear(rad[1],v,0)+(1/r+ct/q)*linear(rad[0],v,0)+linear(ang[1],v,1)/r**2-st*linear(ang[0],v,1)/(r*q)+linear(ang[1],v,2)/q**2
    f,fr,ft,fp,frr,frt,frp,ftt,ftp,fpp=j
    gf=[er[i]*fr+et[i]*ft/r+ep[i]*fp/q for i in range(3)]
    hrr=frr;htt=ftt/r**2+fr/r;hpp=fpp/q**2+ct*fr/q-st*ft/(r*q)
    hrt=frt/r-ft/r**2;hrp=frp/q-ct*fp/q**2;htp=ftp/(r*q)+st*fp/q**2
    hf=[[er[i]*er[k]*hrr+et[i]*et[k]*htt+ep[i]*ep[k]*hpp+(er[i]*et[k]+et[i]*er[k])*hrt+(er[i]*ep[k]+ep[i]*er[k])*hrp+(et[i]*ep[k]+ep[i]*et[k])*htp for k in range(3)] for i in range(3)]
    lf=hrr+htt+hpp
    uc=exact(I.mpf(1)/50);wc=[exact(I.mpf(v)/1000) for v in (10,13,-9)]
    u=1+uc*f;w=[c*f for c in wc]
    vel=[exact(I.mpf(3)/25),exact(-I.mpf(1)/25)];phase=[exact(I.mpf(1)/4),exact(I.mpf(1)/5)]
    rho=exact(I.mpf(189)/200)
    kin=[[1+exact(I.mpf(3)/10)*rho,exact(I.mpf(1)/10)*rho],[exact(I.mpf(1)/10)*rho,1+exact(I.mpf(1)/5)*rho]]
    direction=[et[i]/r-3*ep[i]/q for i in range(3)]
    A=sum(vel[i]*kin[i][k]*vel[k] for i in range(2) for k in range(2))/2+(rho-1)**2
    Bterm=sum(phase[i]*kin[i][k]*phase[k] for i in range(2) for k in range(2))/2*sum(v*v for v in direction)
    S=[-sum(vel[i]*kin[i][k]*phase[k] for i in range(2) for k in range(2))*v for v in direction]
    pi=exact(I.pi)
    def residual(du,dw,div):
        trace=sum(dw[i][i] for i in range(3))
        strain=[[dw[i][k]+dw[k][i]-(exact(I.mpf(2)/3)*trace if i==k else 0) for k in range(3)] for i in range(3)]
        H=-8*u**-5*du-u**-12*sum(v*v for row in strain for v in row)-16*pi*exact(I.mpf(1)/1000)*(A+Bterm*u**-4)
        M=[u**-10*div[i]-8*pi*exact(I.mpf(1)/1000)*u**-4*S[i] for i in range(3)]
        return [H]+M,strain
    dw=[grad(v) for v in w];trace=sum(dw[i][i] for i in range(3))
    strain=[[dw[i][k]+dw[k][i]-(exact(I.mpf(2)/3)*trace if i==k else 0) for k in range(3)] for i in range(3)]
    div=[sum(grad(strain[i][k])[k] for k in range(3)) for i in range(3)]
    discrete,_=residual(lap(u),dw,div)
    edw=[[wc[i]*gf[k] for k in range(3)] for i in range(3)]
    ediv=[wc[i]*lf+sum(wc[k]*hf[i][k] for k in range(3))/3 for i in range(3)]
    continuum,_=residual(uc*lf,edw,ediv)
    e=Evaluator(n,9)
    actual,_=e.residual(u.c,np.array([v.c for v in w]))
    errors=[];bounds=[];direct_bounds=[];discrepancies=[];exact_bounds=[];exact_max=[]
    for i in range(4):
        err=float(np.max(abs(actual[i]-continuum[i].c)))
        # Also covers subtraction and reduction rounding in measuring the error.
        disc=float(np.max(abs(actual[i]-discrete[i].c)))
        db=float(np.max(discrete[i].e));eb=float(np.max(continuum[i].e))
        bound=float(up(disc+db+eb+8*EPS*(np.max(abs(actual[i]))+np.max(abs(continuum[i].c)))))
        errors.append(err);bounds.append(bound);direct_bounds.append(db);discrepancies.append(disc);exact_bounds.append(eb);exact_max.append(float(np.max(abs(continuum[i].c))))
    return dict(n=n,max_abs_error=errors,floating_error_bound=bounds,direct_matrix_radius=direct_bounds,FFT_to_direct_center=discrepancies,continuum_radius=exact_bounds,exact_residual_max=exact_max,error_over_bound=(np.array(errors)/bounds).tolist())

def main():
    rows=[]
    for n in (8,16,32,64):
        row=calculate(n);rows.append(row);print(json.dumps(row),flush=True)
    orders=[]
    for a,b in zip(rows[1:],rows[2:]):
        av=np.array(a['max_abs_error']);bv=np.array(b['max_abs_error']);ua=np.array(a['floating_error_bound']);ub=np.array(b['floating_error_bound'])
        resolved=(av>10*ua)&(bv>10*ub)
        low=np.log2((av-ua)/(bv+ub)) if np.all(av>ua) else np.full(4,np.nan)
        nominal=np.log2(av/bv)
        orders.append(dict(pair=[a['n'],b['n']],resolved_above_10x_bound=resolved.tolist(),nominal=nominal.tolist(),lower_order=[float(x) if np.isfinite(x) else None for x in low],passes_p_1_5=(resolved&(low>=1.5)).tolist()))
    passed=all(all(v['passes_p_1_5']) for v in orders)
    report=dict(status='PASS' if passed else 'FAIL',reference_commit='5799e440b5f27bfe81c3ce5c44814a7d2d652157',components=['H','Mx','My','Mz'],fixture='s=5(r-1/20); F=exp(12(s-1))*(1+sin(theta)/5+3*cos(phi)/20+sin(theta+phi)/10); psi=1+F/50; W=(10,13,-9)*F/1000 (Cartesian)',rows=rows,orders=orders,criterion='Each component on 16->32 and 32->64: both errors >10 times floating bound; conservative lower p >=1.5. No roundoff exemption.',consequence='Proceed to saved-baseline-only audit.' if passed else 'STOP: no baseline residual audit or repair cycle.')
    (OUT/'qualification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('QUALIFICATION',report['status'],flush=True)
    return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
