"""Experimental tensor-divergence solve. Frozen audit is imported unchanged.

The optional angular fourth difference damps collocated angular checkerboards;
its coefficient is O(h_theta**2), so its smooth-field defect is second order.
It does not repair the known first-order composed radial boundary truncation.
"""
import sys, os, json, hashlib, time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator, gmres, spilu
from solve_toroidal_helical_coupled_constraints import Grid, G
from diagnose_toroidal_boundary_convergence import Audit

def output_dir():
    return Path(__file__).parent/'results'/os.environ.get('DEE_EXPERIMENT','consistent-momentum')

class FourierPreconditioner:
    """Exact scalar inversion by periodic phi blocks; no symmetry imposed on W."""
    def __init__(self,L,nphi):
        from scipy.sparse.linalg import splu
        self.nphi=nphi
        blocks=[L[::nphi,k::nphi].tocsc() for k in range(nphi)]
        active=[k for k,b in enumerate(blocks) if b.nnz]
        self.factors=[splu(sum((blocks[k]*np.exp(2j*np.pi*m*k/nphi) for k in active))) for m in range(nphi//2+1)]
    def solve(self,x):
        transformed=np.fft.rfft(x.reshape(-1,self.nphi),axis=1)
        solved=np.column_stack([lu.solve(transformed[:,m]) for m,lu in enumerate(self.factors)])
        return np.fft.irfft(solved,n=self.nphi,axis=1).ravel()

def scalar_flux_laplacian(g):
    """Conservative face fluxes identical in algebra to the frozen H audit.

    Add negative angular fourth differences to suppress grid-scale oscillation.
    Smooth-field correction vanishes quadratically as angular spacing decreases.
    """
    from validate_toroidal_helical_dee_source import R,RMIN,RMAX
    nr,nt,np_=g.shape;hr=(RMAX-RMIN)/(nr-1);ht=2*np.pi/nt;hp=2*np.pi/np_
    ids=np.arange(g.size).reshape(g.shape);rows=[];cols=[];vals=[]
    for i,j,k in np.ndindex(g.shape):
        if i in (0,nr-1):continue
        row=ids[i,j,k];r=g.r[row];t=g.t[row];q=g.q[row]
        terms=[(ids[i+1,j,k],(r+hr/2)*(R+(r+hr/2)*np.cos(t))/(hr**2*r*q)),(ids[i-1,j,k],(r-hr/2)*(R+(r-hr/2)*np.cos(t))/(hr**2*r*q)),(ids[i,(j+1)%nt,k],(R+r*np.cos(t+ht/2))/(ht**2*r*r*q)),(ids[i,(j-1)%nt,k],(R+r*np.cos(t-ht/2))/(ht**2*r*r*q)),(ids[i,j,(k+1)%np_],1/(hp**2*q*q)),(ids[i,j,(k-1)%np_],1/(hp**2*q*q))]
        terms.append((row,-sum(v for _,v in terms)))
        for axis,h,scale in [(1,ht,r*r),(2,hp,q*q)]:
            for offset,weight in [(-2,1),(-1,-4),(0,6),(1,-4),(2,1)]:
                index=[i,j,k];index[axis]=(index[axis]+offset)%g.shape[axis]
                terms.append((ids[tuple(index)],-.25*weight/(h*h*scale)))
        for col,value in terms:rows.append(row);cols.append(col);vals.append(value)
    return sp.csr_matrix((vals,(rows,cols)),shape=(g.size,g.size))

def scalar_fourth_laplacian(g):
    """Fourth-order chart Laplacian, using six points for radial D2 closures."""
    nr,nt,np_=g.shape
    from validate_toroidal_helical_dee_source import RMIN,RMAX
    def matrix(n,h,degree,periodic=False):
        a=sp.lil_matrix((n,n));count=5 if degree==1 or periodic else 6
        for i in range(n):
            js=np.arange(i-2,i+3) if periodic else np.arange(min(max(i-count//2,0),n-count),min(max(i-count//2,0),n-count)+count)
            offsets=(js-i).astype(float);rhs=np.zeros(count);rhs[degree]=1 if degree==1 else 2
            weights=np.linalg.solve(np.array([offsets**k for k in range(count)]),rhs)/h**degree
            for j,w in zip(js,weights):a[i,j%n]+=w
        return a.tocsr()
    ir,it,ip=[sp.eye(n,format='csr') for n in g.shape]
    kron=lambda a,b,c:sp.kron(sp.kron(a,b,format='csr'),c,format='csr')
    hr=(RMAX-RMIN)/(nr-1);ht=2*np.pi/nt;hp=2*np.pi/np_
    dr=kron(matrix(nr,hr,1),it,ip);dt=kron(ir,matrix(nt,ht,1,True),ip)
    return (kron(matrix(nr,hr,2),it,ip)+sp.diags(1/g.r+np.cos(g.t)/g.q)@dr+kron(ir,matrix(nt,ht,2,True),ip).multiply((1/g.r**2)[:,None])-sp.diags(np.sin(g.t)/(g.r*g.q))@dt+sp.diags(1/g.q**2)@kron(ir,it,matrix(np_,hp,2,True))).tocsr()

class TensorGrid(Grid):
    def __init__(self,n,stabilization=0.25):
        super().__init__(n+1,n,n,build_vector=False)
        self.stabilization=stabilization
        self.history=[]
        h=2*np.pi/n
        # Periodic fourth differences: negative definite correction to Delta_L.
        def fourth(v,axis):
            return (np.roll(v,2,axis)-4*np.roll(v,1,axis)+6*v-4*np.roll(v,-1,axis)+np.roll(v,-2,axis))/h**2
        def apply(x):
            w=np.zeros((3,self.size));w[:,self.idx]=x.reshape(3,-1)
            a=self.longitudinal(w)
            y=np.array([sum(self.grad[j]@a[i,j] for j in range(3)) for i in range(3)])
            if stabilization:
                for i in range(3):
                    v=w[i].reshape(self.shape)
                    y[i]-=stabilization*(fourth(v,1).ravel()/self.r**2+fourth(v,2).ravel()/self.q**2)
            return y[:,self.idx].ravel()
        m=3*len(self.idx)
        self.vector=LinearOperator((m,m),matvec=apply,dtype=float)
        if os.environ.get('DEE_SCALAR_ORDER','2')=='4':
            self.lap=scalar_fourth_laplacian(self)
            self.L=(self.lap[self.idx]@self.inject).tocsc()
        if os.environ.get('DEE_SCALAR_ORDER')=='flux':
            self.lap=scalar_flux_laplacian(self)
            self.L=(self.lap[self.idx]@self.inject).tocsc()
        self.ilu=FourierPreconditioner(self.L,n) if os.environ.get('DEE_PRECONDITIONER')=='fourier' else spilu(self.L,drop_tol=1e-4,fill_factor=8)

    def solve(self):
        u=np.ones(self.size)
        for iteration in range(30):
            w=self.solve_vector(u);a=self.longitudinal(w);a2=np.sum(a*a,axis=(0,1))
            for newton in range(12):
                e=self.energy(u);f=self.lap@(u-1)+a2*u**-7/8+2*np.pi*G*e*u**5;res=f[self.idx]
                if np.max(abs(res))<1e-10:break
                de=-2*self.spatial*u**-5
                jac=self.L+sp.diags((-7*a2*u**-8/8+2*np.pi*G*(de*u**5+5*e*u**4))[self.idx])
                pre=LinearOperator(jac.shape,matvec=self.ilu.solve,dtype=float)
                step,info=gmres(jac,-res,M=pre,rtol=1e-11,atol=1e-14,restart=40,maxiter=100)
                if info or np.linalg.norm(jac@step+res)>max(1e-13,1e-10*np.linalg.norm(res)):raise RuntimeError('Newton linear solve failed')
                alpha=1.
                while np.min(u[self.idx]+alpha*step)<=0:alpha/=2
                u[self.idx]+=alpha*step
            else:raise RuntimeError('Newton failed')
            wnew=self.solve_vector(u);change=np.max(abs(wnew-w));w=wnew
            a2=np.sum(self.longitudinal(w)**2,axis=(0,1))
            h=(self.lap@(u-1)+a2*u**-7/8+2*np.pi*G*self.energy(u)*u**5)[self.idx]
            checkpoint=output_dir()/f'checkpoint_{self.shape[1]}_{iteration}.npz'
            np.savez_compressed(checkpoint,psi=u,W_cartesian=w,Hmax=np.max(abs(h)),change=change)
            if change<1e-12 and np.max(abs(h))<1e-9:return u,w,iteration+1
        raise RuntimeError('coupled iteration failed')

    def solve_vector(self,u):
        rhs=(8*np.pi*G*u[self.idx]**6*self.S[:,self.idx]).ravel()
        n=len(self.idx)
        pre=LinearOperator(self.vector.shape,matvec=lambda x:np.concatenate([self.ilu.solve(v) for v in x.reshape(3,n)]),dtype=float)
        hist=[]
        x,info=gmres(self.vector,rhs,M=pre,rtol=1e-10,atol=1e-13,restart=60,maxiter=100,callback=lambda r:hist.append(float(r)),callback_type='pr_norm')
        residual=float(np.linalg.norm(self.vector@x-rhs))
        self.history.append(dict(info=int(info),true_residual=residual,rhs_norm=float(np.linalg.norm(rhs)),history=hist))
        print(json.dumps(dict(gmres_info=int(info),iterations=len(hist),true_residual=residual)),flush=True)
        if info: raise RuntimeError('GMRES failed')
        w=np.zeros((3,self.size));w[:,self.idx]=x.reshape(3,n)
        return w

def main():
    root=Path(__file__).resolve().parents[1];out=output_dir()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Preserving existing evidence: choose a fresh DEE_EXPERIMENT instead of {out}')
    out.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).read_bytes();(out/'executed_solver.py').write_bytes(source)
    (out/'settings.json').write_text(json.dumps(dict(solver_sha256=hashlib.sha256(source).hexdigest(),levels=os.environ.get('DEE_LEVELS','8,16,32,64'),scalar_order=os.environ.get('DEE_SCALAR_ORDER','2'),preconditioner=os.environ.get('DEE_PRECONDITIONER','ilu'),stabilization=.25,gmres_rtol=1e-10,gmres_atol=1e-13,numpy=np.__version__,scipy=__import__('scipy').__version__),indent=2)+'\n')
    audit=root/'reproduce/diagnose_toroidal_boundary_convergence.py'
    assert hashlib.sha1(b'blob '+str(audit.stat().st_size).encode()+b'\0'+audit.read_bytes()).hexdigest()=='5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
    rows=[]
    for n in map(int,os.environ.get('DEE_LEVELS','8,16,32,64').split(',')):
        start=time.time();g=TensorGrid(n);u,w,it=g.solve()
        path=out/f'fields_{n}.npz';np.savez_compressed(path,psi=u.reshape(g.shape),W_cartesian=w.reshape((3,)+g.shape))
        row=dict(resolution=list(g.shape),iterations=it,seconds=time.time()-start,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),audit=Audit(g.shape).summarize(u.reshape(g.shape),w.reshape((3,)+g.shape)),gmres=g.history,diagnostics=g.diagnostics(u,w))
        (out/f'result_{n}.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
        print(json.dumps(dict(completed=n,seconds=row['seconds'])),flush=True)
    checks=[]
    for band in ('inner_boundary_band','outer_boundary_band'):
        for metric in ('independent_flux_H','independent_tensor_M_order2'):
            for a,b in zip(rows,rows[1:]):
                order=float(np.log2(a['audit'][band][metric]['volume_weighted_rms']/b['audit'][band][metric]['volume_weighted_rms']))
                checks.append(dict(band=band,metric=metric,pair=[a['resolution'][1],b['resolution'][1]],order=order,passed=bool(np.isfinite(order) and order>=1.5)))
    report=dict(status='PASS' if len(checks)==12 and all(c['passed'] for c in checks) else 'FAIL',checks=checks,passed=sum(c['passed'] for c in checks),rows=rows)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2),flush=True)
if __name__=='__main__':main()
