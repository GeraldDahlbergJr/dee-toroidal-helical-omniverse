"""Boundary audit of frozen Run42 fields; no equations, seed, or solves changed.

Separate algebraic solve residuals from physical tensor-divergence and
face-flux residuals. Evaluate fixed physical radial bands, including an
independent fourth-order derivative audit at the endpoints.
"""
from __future__ import annotations
import hashlib
import json
import zipfile
from pathlib import Path
import numpy as np
from dee_stress_energy_3p1 import coeffs, potential
from validate_toroidal_helical_dee_source import R,RMIN,RMAX,M,RHO,PI,KTH,KPS

G=1e-3

def radial_matrix(n,h,order=2):
    count=3 if order==2 else 5
    out=np.zeros((n,n))
    for i in range(n):
        start=min(max(i-count//2,0),n-count)
        js=np.arange(start,start+count); offsets=js-i
        powers=np.array([offsets.astype(float)**k for k in range(count)])
        rhs=np.zeros(count); rhs[1]=1.
        out[i,js]=np.linalg.solve(powers,rhs)/h
    return out

class Audit:
    def __init__(self,shape):
        nr,nt,nph=shape
        self.shape=tuple(shape); self.dr=(RMAX-RMIN)/(nr-1);self.dt=2*np.pi/nt;self.dp=2*np.pi/nph
        self.r=np.linspace(RMIN,RMAX,nr)[:,None,None]
        self.t=np.arange(nt)[None,:,None]*self.dt; self.p=np.arange(nph)[None,None,:]*self.dp
        self.q=R+self.r*np.cos(self.t)
        broadcast=lambda *v:np.array(np.broadcast_arrays(*v))
        self.er=broadcast(np.cos(self.t)*np.cos(self.p)+0*self.r,np.cos(self.t)*np.sin(self.p)+0*self.r,np.sin(self.t)+0*self.p+0*self.r)
        self.et=broadcast(-np.sin(self.t)*np.cos(self.p)+0*self.r,-np.sin(self.t)*np.sin(self.p)+0*self.r,np.cos(self.t)+0*self.p+0*self.r)
        self.ep=broadcast(-np.sin(self.p)+0*self.r+0*self.t,np.cos(self.p)+0*self.r+0*self.t,0*self.r+0*self.t+0*self.p)
        self.drad={order:radial_matrix(nr,self.dr,order) for order in (2,4)}
        kt,kp,la=coeffs(RHO); kin=np.array([[1.,0,0],[0,kt,la],[0,la,kp]])
        directions=self.et/self.r-M*self.ep/self.q
        D=np.zeros((3,3)+self.shape);D[1]=KTH*directions;D[2]=KPS*directions
        self.spatial=np.einsum('ab,ai...,bi...->...',kin,D,D)
        self.temporal=float(PI@kin@PI);self.S=-np.einsum('ab,a,bi...->i...',kin,PI,D)

    def angular_derivative(self,u,axis,h,order):
        if order==2:return (np.roll(u,-1,axis)-np.roll(u,1,axis))/(2*h)
        return (-np.roll(u,-2,axis)+8*np.roll(u,-1,axis)-8*np.roll(u,1,axis)+np.roll(u,2,axis))/(12*h)

    def grad(self,u,order):
        ur=np.einsum('ij,jkl->ikl',self.drad[order],u)
        ut=self.angular_derivative(u,1,self.dt,order);up=self.angular_derivative(u,2,self.dp,order)
        return self.er*ur+self.et*ut/self.r+self.ep*up/self.q

    def tensor(self,w,order):
        dw=np.array([self.grad(v,order) for v in w]);div=np.einsum('ii...->...',dw)
        return dw+dw.swapaxes(0,1)-2/3*np.eye(3)[:,:,None,None,None]*div

    def momentum(self,u,w,order):
        a=self.tensor(w,order)
        da=np.array([sum(self.grad(a[i,j],order)[j] for j in range(3)) for i in range(3)])
        return u**-10*da-8*np.pi*G*u**-4*self.S

    def flux_hamiltonian(self,u,w):
        r,t,q=self.r,self.t,self.q;dr,dt,dp=self.dr,self.dt,self.dp
        rf=(r[1:]+r[:-1])/2;fr=rf*(R+rf*np.cos(t))*(u[1:]-u[:-1])/dr
        lap=(fr[1:]-fr[:-1])/dr/(r[1:-1]*q[1:-1])
        ft=(R+r*np.cos(t+dt/2))/r*(np.roll(u,-1,1)-u)/dt
        lap+=((ft-np.roll(ft,1,1))/dt/(r*q))[1:-1]
        lap+=((np.roll(u,-1,2)-2*u+np.roll(u,1,2))/(dp*dp*q*q))[1:-1]
        a=self.tensor(w,2);a2=np.einsum('ij...,ij...->...',a,a)
        en=.5*self.temporal+potential(RHO)+.5*self.spatial*u**-4
        result=np.full(self.shape,np.nan)
        result[1:-1]=-8*u[1:-1]**-5*lap-a2[1:-1]*u[1:-1]**-12-16*np.pi*G*en[1:-1]
        return result

    def norms(self,v,u,mask):
        scalar=v.ndim==3
        sq=v*v if scalar else np.sum(v*v,axis=0)
        weight=np.broadcast_to(self.r*self.q,self.shape)*u**6
        if not np.any(mask):raise ValueError('empty diagnostic region')
        return {'volume_weighted_rms':float(np.sqrt(np.sum((sq*weight)[mask])/np.sum(weight[mask]))),
                'max_abs' if scalar else 'max_vector_norm':float(np.sqrt(np.max(sq[mask])))}

    def summarize(self,u,w):
        rv=self.r[:,0,0];interior=np.ones(self.shape,dtype=bool);interior[[0,-1]]=False
        bands={'whole_interior':[RMIN,RMAX],'inner_boundary_band':[RMIN,.10],'bulk':[.10,.20],'outer_boundary_band':[.20,RMAX]}
        mom={k:self.momentum(u,w,k) for k in (2,4)}; H=self.flux_hamiltonian(u,w)
        out={}
        for name,(lo,hi) in bands.items():
            radial=(rv>=lo-1e-12)&(rv<=hi+1e-12)
            mask=np.broadcast_to(radial[:,None,None],self.shape)&interior
            out[name]={'r_interval':[lo,hi],'interior_radial_samples':int(np.sum(radial[1:-1])),
                'independent_flux_H':self.norms(H,u,mask),
                'independent_tensor_M_order2':self.norms(mom[2],u,mask),
                'independent_tensor_M_order4':self.norms(mom[4],u,mask)}
        for index,name in [(0,'inner_endpoint'),(-1,'outer_endpoint')]:
            mask=np.zeros(self.shape,dtype=bool);mask[index]=True
            out[name]={'r':float(rv[index]),'independent_tensor_M_order2':self.norms(mom[2],u,mask),
                'independent_tensor_M_order4':self.norms(mom[4],u,mask),
                'note':'one-sided PDE compatibility audit; Dirichlet boundary values themselves do not enforce the constraints at endpoints'}
        return out

def main():
    root=Path(__file__).parent;frozen=root/'frozen/run42';manifest=json.loads((frozen/'manifest.json').read_text())
    archive=frozen/'artifacts.zip'
    assert 'sha256:'+hashlib.sha256(archive.read_bytes()).hexdigest()==manifest['source']['artifact_digest']
    rows=[]
    with zipfile.ZipFile(archive) as z:
        checkpoint=json.loads(z.read('toroidal_helical_coupled_constraints.json'))
        for old in checkpoint['rows']:
            nr,nt,nph=old['resolution'];name=f'toroidal_coupled_fields_{nr}_{nt}.npz'
            raw=z.read(name);assert hashlib.sha256(raw).hexdigest()==checkpoint['raw_field_sha256'][name]
            import io
            with np.load(io.BytesIO(raw)) as data:u,w=data['psi'],data['W_cartesian']
            assert np.all(np.isfinite(u)) and np.all(np.isfinite(w)) and np.min(u)>0
            audit=Audit(u.shape);regions=audit.summarize(u,w)
            # Reproduce the old exactly-defined bulk norm before adding new bands.
            rv=audit.r[:,0,0]; mask=np.broadcast_to(((rv>=.10)&(rv<=.20))[:,None,None],u.shape)
            oldnorm=audit.norms(audit.momentum(u,w,2),u,mask)['volume_weighted_rms']
            assert np.isclose(oldnorm,old['independent_M_vector_bulk_rms'],rtol=1e-9,atol=1e-13)
            rows.append({'resolution':[nr,nt,nph],'field_sha256':hashlib.sha256(raw).hexdigest(),
                'archived_bulk_norm_reproduced':True,'regions':regions})
            print(f'Audited {nr} x {nt} x {nph}',flush=True)
    orders={}
    for region in rows[0]['regions']:
        orders[region]={}
        for metric in rows[0]['regions'][region]:
            if not metric.startswith('independent_'):continue
            values=[r['regions'][region][metric]['volume_weighted_rms'] for r in rows]
            orders[region][metric]=[float(np.log2(values[i]/values[i+1])) for i in range(2)]
    checks={region:all(p>=1.5 for metric in ('independent_flux_H','independent_tensor_M_order2') for p in orders[region][metric]) for region in ('inner_boundary_band','outer_boundary_band')}
    out={'checkpoint_name':'Frozen Run42 boundary-convergence audit','status':'BOUNDARY_GATE_PASS' if all(checks.values()) else 'BOUNDARY_GATE_NOT_MET',
        'diagnostic_integrity':'PASS','source_run':42,'frozen_source_commit':manifest['source']['source_commit'],
        'frozen_archive_digest':manifest['source']['artifact_digest'],'diagnostic_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'preserved':'Run42 equations, matter seed, G, geometry, boundary conditions, source solver and raw fields; no new solve or parameter change',
        'method':'independent physical tensor divergence with second- and fourth-order first derivatives; conservative face-flux Hamiltonian; fixed physical radial bands; volume weights sqrt(gamma)',
        'boundary_criterion':'both independent second-order H and M RMS must show order >=1.5 on both refinement pairs in each boundary band',
        'rows':rows,'observed_orders':orders,'boundary_gate':checks,
        'limitations':['Three frozen resolutions only; endpoint derivatives are extrapolation-based compatibility diagnostics.','Fourth-order audit does not upgrade second-order solved fields or substitute for the primary second-order convergence gate.','Axisymmetric stress seed remains unchanged; no generic nonaxisymmetric or evolved solution claimed.']}
    path=root/'toroidal_boundary_convergence.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'boundary_gate':checks,'orders':orders},indent=2))
    # A valid diagnostic can discover failure. Do not turn scientific gate failure
    # into a missing artifact or pretend that a successful workflow validates it.
if __name__=='__main__':main()
