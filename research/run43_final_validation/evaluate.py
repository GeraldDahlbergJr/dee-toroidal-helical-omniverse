"""Standalone saved-array continuum audit. No project imports, solves or writes to inputs."""
from pathlib import Path
import hashlib, json, math
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'research/results/coordinated-replay-20261103-1831/baseline'

def radial(n,width,degree):
    # Independent local polynomial interpolation, differentiated at each node.
    mat=np.zeros((n+1,n+1))
    # Exact uniform-grid offsets avoid cancellation in floating polynomial
    # normalization. Convert only the final physical derivative weights.
    h=sp.Rational(1,5*n)
    for i in range(n+1):
        start=min(max(i-width//2,0),n+1-width)
        ids=list(range(start,start+width))
        weights=sp.finite_diff_weights(degree,[k-i for k in ids],0)[degree][-1]
        for k,weight in zip(ids,weights):
            mat[i,k]=float(weight/h**degree)
    return mat

class Evaluator:
    def __init__(self,n,width,cutoff=None):
        self.cutoff=cutoff
        self.n=n; self.shape=(n+1,n,n)
        self.r=np.linspace(.05,.25,n+1)[:,None,None]
        self.t=2*np.pi*np.arange(n)[None,:,None]/n
        self.p=2*np.pi*np.arange(n)[None,None,:]/n
        self.q=1+self.r*np.cos(self.t)
        self.frame=np.array(np.broadcast_arrays(np.cos(self.t)*np.cos(self.p)+0*self.r,np.cos(self.t)*np.sin(self.p)+0*self.r,np.sin(self.t)+0*self.p+0*self.r))
        self.et=np.array(np.broadcast_arrays(-np.sin(self.t)*np.cos(self.p)+0*self.r,-np.sin(self.t)*np.sin(self.p)+0*self.r,np.cos(self.t)+0*self.p+0*self.r))
        self.ep=np.array(np.broadcast_arrays(-np.sin(self.p)+0*self.r+0*self.t,np.cos(self.p)+0*self.r+0*self.t,0*self.r+0*self.t+0*self.p))
        self.d1=radial(n,width,1);self.d2=radial(n,width,2)
        self.xyz=np.array(np.broadcast_arrays(self.q*np.cos(self.p),self.q*np.sin(self.p),self.r*np.sin(self.t)))
        # Independently reconstruct kinetic matrix and normal/spatial projections.
        kin=np.array([[1+.30*.945,.10*.945],[.10*.945,1+.20*.945]])
        velocity=np.array([.12,-.04]);phase=np.array([.25,.20])
        direction=self.et/self.r-3*self.ep/self.q
        self.A=.5*velocity@kin@velocity+(.945-1)**2
        self.B=.5*(phase@kin@phase)*np.sum(direction**2,axis=0)
        self.S=-(velocity@kin@phase)*direction
    def angular(self,u,axis,degree=1):
        k=np.fft.fftfreq(self.n,1/self.n); shape=[1]*3;shape[axis]=self.n
        factor=(1j*k)**degree
        if self.cutoff is not None: factor[abs(k)>self.cutoff]=0
        if degree%2: factor[self.n//2]=0 # real periodic interpolant Nyquist derivative at nodes
        return np.fft.ifft(np.fft.fft(u,axis=axis)*factor.reshape(shape),axis=axis).real
    def grad(self,u):
        return self.frame*np.einsum('ij,jkl->ikl',self.d1,u)+self.et*self.angular(u,1)/self.r+self.ep*self.angular(u,2)/self.q
    def lap(self,u):
        return (np.einsum('ij,jkl->ikl',self.d2,u)+(1/self.r+np.cos(self.t)/self.q)*np.einsum('ij,jkl->ikl',self.d1,u)+self.angular(u,1,2)/self.r**2-np.sin(self.t)*self.angular(u,1)/(self.r*self.q)+self.angular(u,2,2)/self.q**2)
    def residual(self,u,w):
        dw=np.array([self.grad(v) for v in w]); trace=np.einsum('ii...->...',dw)
        strain=dw+dw.swapaxes(0,1)-2/3*np.eye(3)[:,:,None,None,None]*trace
        div=np.array([sum(self.grad(strain[i,j])[j] for j in range(3)) for i in range(3)])
        hterms=np.array([-8*u**-5*self.lap(u),-u**-12*np.sum(strain**2,axis=(0,1)),-16*np.pi*.001*(self.A+self.B*u**-4)])
        mterms=np.array([u**-10*div,-8*np.pi*.001*u**-4*self.S])
        return np.concatenate([hterms.sum(axis=0)[None],mterms.sum(axis=0)]), np.concatenate([np.sum(abs(hterms),axis=0)[None],np.sum(abs(mterms),axis=0)])
    def regions(self):
        r=self.r[:,0,0];return {'interior':np.arange(1,self.n),'inner_band':np.where((r>.05)&(r<=.10000000001))[0],'outer_band':np.where((r>=.19999999999)&(r<.25))[0],'inner_first_row':np.array([1]),'outer_first_row':np.array([self.n-1])}
    def rms(self,v,u,ids):
        wt=np.broadcast_to(self.r*self.q,self.shape)*u**6
        return np.sqrt(np.sum(v[:,ids]**2*wt[ids],axis=(1,2,3))/np.sum(wt[ids]))

def fixture():
    # One analytic Cartesian polynomial field; toroidal pullback has radial,
    # both angular, and mixed derivatives. Exact derivatives use Cartesian symbols.
    x,y,z=sp.symbols('x y z');coords=(x,y,z)
    u=1+sp.Rational(1,1000)*(x*x*y+y*z*z+z*x)
    w=sp.Matrix([x*x*y+y*z*z,x*y*y+x*z*z,x*y*z+z**3])/1000
    jac=w.jacobian(coords); strain=jac+jac.T-sp.eye(3)*sp.trace(jac)*sp.Rational(2,3)
    lap=sum(sp.diff(u,c,2) for c in coords)
    div=sp.Matrix([sum(sp.diff(strain[i,j],coords[j]) for j in range(3)) for i in range(3)])
    def sample(expr,e): return np.broadcast_to(np.asarray(sp.lambdify(coords,expr,'numpy')(*e.xyz),float),e.shape)
    rows=[]
    for n in (8,16,32,64):
        e=Evaluator(n,min(9,n+1)); us=sample(u,e);ws=np.array([sample(v,e) for v in w])
        a2=sample(sum(v*v for v in strain),e)
        exact=np.concatenate([(-8*us**-5*sample(lap,e)-us**-12*a2-16*np.pi*.001*(e.A+e.B*us**-4))[None],us**-10*np.array([sample(v,e) for v in div])-8*np.pi*.001*us**-4*e.S])
        actual,_=e.residual(us,ws)
        rows.append({'n':n,'max_abs_error':np.max(abs(actual-exact),axis=(1,2,3)).tolist(),'exact_residual_max':np.max(abs(exact),axis=(1,2,3)).tolist()})
    return {'psi':str(u),'W':[str(v) for v in w],'method':'Exact symbolic Cartesian derivatives, compared with toroidal polynomial/Fourier derivatives. No manufactured solve. Errors reported without fitted tolerance.','rows':rows}

def main():
    paths=[ROOT/'research/RUN43_CONSOLIDATED_MATHEMATICAL_SPECIFICATION_20261106.md',ROOT/'docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md',ROOT/'reproduce/dee_stress_energy_3p1.py',ROOT/'reproduce/validate_toroidal_helical_dee_source.py',ROOT/'reproduce/solve_toroidal_helical_coupled_constraints.py',ROOT/'reproduce/diagnose_toroidal_boundary_convergence.py',ROOT/'research/audit_run43_arrays.py',ROOT/'research/RUN43_COORDINATED_REPLAY_20261103.md']
    paths+=list(BASE.glob('*')); paths+=[ROOT/'research/results/coordinated-replay-20261103-1831/coordination/frozen_audit.py',ROOT/'research/results/coordinated-replay-20261103-1831/validation/baseline-independent-validation.json']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
    (OUT/'input_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    result={'qualification_status':'FAILED: see qualification.json; diagnostic values are unqualified','components':['H','Mx','My','Mz'],'fixture':fixture(),'rows':[],'orders':[]}
    for n in (8,16,32,64):
        path=BASE/f'fields_{n}.npz';expected=json.loads((BASE/f'result_{n}.json').read_text())['sha256'];assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
        with np.load(path) as data:u=data['psi'];w=data['W_cartesian']
        assert u.shape==(n+1,n,n) and w.shape==(3,n+1,n,n)
        assert np.isfinite(u).all() and np.isfinite(w).all() and u.min()>0
        e=Evaluator(n,9); high,scale=e.residual(u,w)
        mid,_=Evaluator(n,7).residual(u,w);low,_=Evaluator(n,5).residual(u,w)
        angular,_=Evaluator(n,9,n//3).residual(u,w)
        uncertainty=np.maximum.reduce([abs(high-mid),abs(mid-low),abs(high-angular)])
        boundary={'psi_trace_max_error':float(np.max(abs(u[[0,-1]]-1))),'W_trace_max_error':float(np.max(abs(w[:,[0,-1]]))),'angular_seams':'Endpoint-exclusive arrays: no independent duplicate 2pi trace is saved; equality of physical seam traces NOT YET ESTABLISHED. Fourier continuation is periodic by construction, not evidence of saved endpoint equality.','adjacent_wrap_max':{}}
        for axis in (1,2):
            boundary['adjacent_wrap_max'][str(axis)]={'psi':float(np.max(abs(np.take(u,0,axis)-np.take(u,-1,axis)))),'W':float(np.max(abs(np.take(w,0,axis+1)-np.take(w,-1,axis+1))))}
        row={'n':n,'boundary':boundary,'regions':{}}
        for name,ids in e.regions().items():
            v=e.rms(high,u,ids);err=e.rms(uncertainty,u,ids);den=e.rms(scale,u,ids)
            row['regions'][name]={'rms':v.tolist(),'max_abs':np.max(abs(high[:,ids]),axis=(1,2,3)).tolist(),'term_scale_rms':den.tolist(),'normalized_rms':(v/den).tolist(),'uncertainty_rms':err.tolist(),'uncertainty_fraction':(err/v).tolist(),'within_10_percent':(err<=.1*v).tolist()}
        result['rows'].append(row)
        print('evaluated',n,flush=True)
    for name in result['rows'][0]['regions']:
        for a,b in zip(result['rows'],result['rows'][1:]):
            av=a['regions'][name];bv=b['regions'][name]
            aa=np.array(av['rms']);bb=np.array(bv['rms']);ua=np.array(av['uncertainty_rms']);ub=np.array(bv['uncertainty_rms'])
            p=np.log2(aa/bb)
            bound=np.full(4,-np.inf); valid=aa>ua
            bound[valid]=np.log2((aa[valid]-ua[valid])/(bb[valid]+ub[valid]))
            result['orders'].append({'region':name,'pair':[a['n'],b['n']],'nominal':p.tolist(),'conservative':[float(v) if np.isfinite(v) else None for v in bound],'passes_p_1_5':(bound>=1.5).tolist(),'both_uncertainties_within_budget':(np.maximum(ua/aa,ub/bb)<=.1).tolist()})
    result['uncertainty_scope']='Empirical envelope of 5/7/9-point radial reconstructions and removal of angular modes above n//3. Null conservative order means uncertainty contains zero; no finite lower order established. Not a rigorous error bound; unresolved angular aliases cannot be bounded from saved samples alone.'
    result['absolute_accuracy']='NOT YET ESTABLISHED: no predeclared absolute physical residual tolerance recovered; algebraic solver tolerances are not physical tolerances.'
    (OUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
