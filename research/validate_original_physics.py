"""Supplementary original-source audit; never replaces the frozen twelve checks.

Only solved arrays are evaluated. Radial seven-point polynomial derivatives and
Fourier angular derivatives differ from both solver and frozen diagnostic.
No synthetic solutions, source adjustments, or new acceptance thresholds.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reproduce'))
from diagnose_toroidal_boundary_convergence import Audit
from dee_stress_energy_3p1 import direct_tensor_check

def radial(n,h,degree):
    out=np.zeros((n,n));m=min(7,n)
    for i in range(n):
        start=min(max(i-m//2,0),n-m);js=np.arange(start,start+m)
        v=js-i;rhs=np.zeros(m);rhs[degree]=1 if degree==1 else 2
        out[i,js]=np.linalg.solve(np.array([v.astype(float)**k for k in range(m)]),rhs)/h**degree
    return out

class Independent:
    def __init__(self,n):
        self.n=n;self.r=np.linspace(.05,.25,n+1)[:,None,None]
        self.t=np.arange(n)[None,:,None]*2*np.pi/n
        self.p=np.arange(n)[None,None,:]*2*np.pi/n
        self.q=1+self.r*np.cos(self.t)
        self.d1=radial(n+1,.2/n,1);self.d2=radial(n+1,.2/n,2)
        bc=lambda *v:np.array(np.broadcast_arrays(*v))
        self.er=bc(np.cos(self.t)*np.cos(self.p)+0*self.r,np.cos(self.t)*np.sin(self.p)+0*self.r,np.sin(self.t)+0*self.p+0*self.r)
        self.et=bc(-np.sin(self.t)*np.cos(self.p)+0*self.r,-np.sin(self.t)*np.sin(self.p)+0*self.r,np.cos(self.t)+0*self.p+0*self.r)
        self.ep=bc(-np.sin(self.p)+0*self.r+0*self.t,np.cos(self.p)+0*self.r+0*self.t,0*self.r+0*self.p+0*self.t)
    def angle(self,v,axis,degree=1):
        k=np.fft.fftfreq(self.n,1/self.n)
        if degree==1:k[self.n//2]=0
        shape=[1]*3;shape[axis]=self.n
        return np.fft.ifft(np.fft.fft(v,axis=axis)*(1j*k.reshape(shape))**degree,axis=axis).real
    def dr(self,v):return np.einsum('ij,jkl->ikl',self.d1,v)
    def grad(self,v):return self.er*self.dr(v)+self.et*self.angle(v,1)/self.r+self.ep*self.angle(v,2)/self.q
    def residuals(self,u,w):
        dw=np.array([self.grad(v) for v in w]);div=np.trace(dw,axis1=0,axis2=1)
        a=dw+dw.swapaxes(0,1)-2/3*np.eye(3)[:,:,None,None,None]*div
        da=np.array([sum(self.grad(a[i,j])[j] for j in range(3)) for i in range(3)])
        lap=np.einsum('ij,jkl->ikl',self.d2,u-1)+(1/self.r+np.cos(self.t)/self.q)*self.dr(u-1)+self.angle(u-1,1,2)/self.r**2-np.sin(self.t)/(self.r*self.q)*self.angle(u-1,1)+self.angle(u-1,2,2)/self.q**2
        # Reconstruct the original action's 4D stress tensor on actual metrics.
        # Vectorize T00 and -Ti0 directly, retaining the mixed kinetic term.
        kin=np.array([[1,0,0],[0,1+.3*.945,.1*.945],[0,.1*.945,1+.2*.945]])
        pi=np.array([0,.12,-.04]);D=np.zeros((3,3)+u.shape)
        direction=self.et/self.r-3*self.ep/self.q;D[1]=.25*direction;D[2]=.20*direction
        spatial=np.einsum('ab,ai...,bi...->...',kin,D,D)*u**-4
        kinetic=-float(pi@kin@pi)+spatial;lm=-.5*kinetic-(.945-1)**2
        energy=float(pi@kin@pi)-lm
        S=-np.einsum('ab,a,bi...->i...',kin,pi,D)
        H=-8*u**-5*lap-u**-12*np.sum(a*a,axis=(0,1))-16*np.pi*.001*energy
        M=u**-10*da-8*np.pi*.001*u**-4*S
        reference=Audit(u.shape)
        expected=.5*reference.temporal+(.945-1)**2+.5*reference.spatial*u**-4
        source_error={'energy_max_abs':float(np.max(abs(energy-expected))),'momentum_max_abs':float(np.max(abs(S-reference.S)))}
        # Compare independent existing 4D construction at fixed real grid samples.
        worst=0.
        for i,j,k in [(0,0,0),(self.n//2,self.n//4,self.n//4),(self.n,self.n//2,self.n//2)]:
            gamma=np.eye(3)*u[i,j,k]**4
            e,s,*_=direct_tensor_check(np.sqrt(.945),pi,D[:,:,i,j,k],gamma)
            worst=max(worst,abs(e-energy[i,j,k]),float(np.max(abs(s-S[:,i,j,k]))))
        source_error['direct_4D_sample_max_abs']=worst
        regions={}
        for name,lo,hi in [('interior',.05,.25),('inner',.05,.10),('bulk',.10,.20),('outer',.20,.25)]:
            mask=np.broadcast_to((self.r>=lo-1e-12)&(self.r<=hi+1e-12),u.shape).copy();mask[[0,-1]]=False
            regions[name]={'H':reference.norms(H,u,mask),'M':reference.norms(M,u,mask)}
        for i,name in [(0,'inner_endpoint'),(-1,'outer_endpoint')]:
            mask=np.zeros(u.shape,bool);mask[i]=True
            regions[name]={'H':reference.norms(H,u,mask),'M':reference.norms(M,u,mask)}
        return dict(n=self.n,source_consistency=source_error,regions=regions)

def main():
    directory=Path(sys.argv[1]);rows=[];arrays={}
    for n in (8,16,32,64):
        p=directory/f'fields_{n}.npz'
        with np.load(p) as z:u,w=z['psi'],z['W_cartesian']
        arrays[n]=(u.copy(),w.copy());row=Independent(n).residuals(u,w)
        row['sha256']=hashlib.sha256(p.read_bytes()).hexdigest();rows.append(row)
        print(json.dumps(row),flush=True)
    differences=[]
    for a,b in [(8,16),(16,32),(32,64)]:
        u,w=arrays[a];v,z=arrays[b];du=u-v[::2,::2,::2];dw=w-z[:,::2,::2,::2]
        audit=Audit(u.shape);mask=np.ones(u.shape,bool);mask[[0,-1]]=False
        differences.append(dict(pair=[a,b],psi=audit.norms(du,u,mask),W=audit.norms(dw,u,mask),psi_max_all=float(np.max(abs(du))),W_max_all=float(np.max(abs(dw)))))
    orders=[]
    for region in rows[0]['regions']:
        for field in ('H','M'):
            vals=[r['regions'][region][field]['volume_weighted_rms'] for r in rows]
            orders.append(dict(region=region,field=field,values=vals,orders=[float(np.log2(a/b)) for a,b in zip(vals,vals[1:])]))
    solution_orders={f:[float(np.log2(a[f]['volume_weighted_rms']/b[f]['volume_weighted_rms'])) for a,b in zip(differences,differences[1:])] for f in ('psi','W')}
    report=dict(scope='Fixed-source conformal initial-data constraints only; supplementary diagnostics with no new pass threshold. Full Einstein-DEE evolution and matter equations not solved by this candidate.',rows=rows,residual_orders=orders,nested_solution_differences=differences,nested_solution_orders=solution_orders)
    (directory/'original-physics-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'solution_orders':solution_orders,'residual_orders':orders},indent=2))
if __name__=='__main__':main()
