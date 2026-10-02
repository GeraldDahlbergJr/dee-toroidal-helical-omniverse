"""Coupled conformal Hamiltonian + vector momentum benchmark on a toroidal shell.

Cartesian components on a toroidal grid avoid coordinate-component connection
bookkeeping. Physical matter fields and normal velocities are held fixed;
energy is recomputed with gamma^-1=psi^-4 delta at every nonlinear iteration.
Run-38 has axisymmetric stress despite its helical phase: this is a 3D grid
benchmark, NOT evidence of a general nonaxisymmetric solve or evolution.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import splu, gmres, LinearOperator
from dee_stress_energy_3p1 import coeffs, potential
from validate_toroidal_helical_dee_source import R, RMIN, RMAX, M, RHO, PI, KTH, KPS

G = 1e-3

def derivative(n, h, periodic=False, second=False):
    a=sp.lil_matrix((n,n))
    for i in range(n):
        if periodic:
            for j,v in ((i-1,1 if second else -1),(i, -2 if second else 0),(i+1,1)):
                a[i,j%n]+=v/(h*h if second else 2*h)
        elif 0<i<n-1:
            for j,v in ((i-1,1 if second else -1),(i,-2 if second else 0),(i+1,1)):
                a[i,j]=v/(h*h if second else 2*h)
        elif not second:
            js=[0,1,2] if i==0 else [n-3,n-2,n-1]
            vs=[-3,4,-1] if i==0 else [1,-4,3]
            for j,v in zip(js,vs): a[i,j]=v/(2*h)
    return a.tocsr()

class Grid:
    def __init__(self,nr,nt,nph,build_vector=True):
        self.shape=(nr,nt,nph)
        r=np.linspace(RMIN,RMAX,nr); t=np.arange(nt)*2*np.pi/nt; p=np.arange(nph)*2*np.pi/nph
        rr,tt,pp=np.meshgrid(r,t,p,indexing='ij'); q=R+rr*np.cos(tt)
        self.r,self.t,self.p,self.q=[x.ravel() for x in (rr,tt,pp,q)]
        er=np.array([np.cos(tt)*np.cos(pp),np.cos(tt)*np.sin(pp),np.sin(tt)]).reshape(3,-1)
        et=np.array([-np.sin(tt)*np.cos(pp),-np.sin(tt)*np.sin(pp),np.cos(tt)]).reshape(3,-1)
        ep=np.array([-np.sin(pp),np.cos(pp),np.zeros_like(pp)]).reshape(3,-1)
        self.basis=(er,et,ep)
        self.size=rr.size
        self.idx=np.arange(self.size).reshape(self.shape)[1:-1].ravel()
        self.inject=sp.csr_matrix((np.ones(len(self.idx)),(self.idx,np.arange(len(self.idx)))),shape=(self.size,len(self.idx)))
        kron=lambda a,b,c:sp.kron(sp.kron(a,b,format='csr'),c,format='csr')
        ir,it,ip=[sp.eye(n,format='csr') for n in self.shape]
        dr=r[1]-r[0]; dt=2*np.pi/nt; dp=2*np.pi/nph
        self.d=[kron(derivative(nr,dr),it,ip),kron(ir,derivative(nt,dt,True),ip),kron(ir,it,derivative(nph,dp,True))]
        d2=[kron(derivative(nr,dr,second=True),it,ip),kron(ir,derivative(nt,dt,True,True),ip),kron(ir,it,derivative(nph,dp,True,True))]
        diag=sp.diags
        self.lap=(d2[0]+diag(1/self.r+np.cos(self.t)/self.q)@self.d[0]+diag(1/self.r**2)@d2[1]-diag(np.sin(self.t)/(self.r*self.q))@self.d[1]+diag(1/self.q**2)@d2[2]).tocsr()
        self.grad=[(diag(er[i])@self.d[0]+diag(et[i]/self.r)@self.d[1]+diag(ep[i]/self.q)@self.d[2]).tocsr() for i in range(3)]
        L=self.lap[self.idx]@self.inject
        self.L=L.tocsc(); self.lu=splu(self.L)
        # Delta_L W = Delta W + 1/3 grad(div W) in Cartesian components.
        # The frozen/default solver still assembles the original block operator.
        # Experimental subclasses that replace self.vector may skip only this
        # otherwise-discarded assembly; all Grid geometry/operators remain identical.
        if build_vector:
            self.vector=sp.bmat([[(L if i==j else sp.csr_matrix(L.shape))+(self.grad[i]@self.grad[j])[self.idx]@self.inject/3 for j in range(3)] for i in range(3)],format='csr')
        k1,k2,c=coeffs(RHO); kinetic=np.array([[1,0,0],[0,k1,c],[0,c,k2]])
        # Chart covectors -> Cartesian covectors using the reciprocal basis.
        D=np.zeros((3,3,self.size)); direction=et/self.r-M*ep/self.q
        D[1]=KTH*direction; D[2]=KPS*direction
        self.temporal=float(PI@kinetic@PI)
        self.spatial=np.einsum('ab,ain,bin->n',kinetic,D,D)
        self.S=-np.einsum('ab,a,bin->in',kinetic,PI,D)
        self.U=potential(RHO)

    def energy(self,u): return .5*self.temporal+.5*self.spatial*u**-4+self.U

    def longitudinal(self,w):
        dw=np.array([[self.grad[i]@w[j] for j in range(3)] for i in range(3)])
        div=sum(dw[i,i] for i in range(3))
        a=dw+dw.swapaxes(0,1)
        for i in range(3): a[i,i]-=2*div/3
        return a

    def solve_vector(self,u):
        rhs=(8*np.pi*G*u[self.idx]**6*self.S[:,self.idx]).ravel()
        n=len(self.idx)
        pre=LinearOperator((3*n,3*n),matvec=lambda x:np.concatenate([self.lu.solve(y) for y in x.reshape(3,n)]))
        x,info=gmres(self.vector,rhs,M=pre,rtol=1e-10,atol=1e-13,restart=60,maxiter=100)
        if info: raise RuntimeError(f'vector GMRES failed: {info}')
        w=np.zeros((3,self.size)); w[:,self.idx]=x.reshape(3,n)
        return w

    def solve(self):
        u=np.ones(self.size)
        for iteration in range(30):
            w=self.solve_vector(u); a=self.longitudinal(w); a2=np.sum(a*a,axis=(0,1))
            # Newton Hamiltonian solve, exact scalar Laplacian Jacobian.
            for _ in range(12):
                e=self.energy(u)
                f=self.lap@(u-1)+a2*u**-7/8+2*np.pi*G*e*u**5
                residual=f[self.idx]
                if np.max(abs(residual))<1e-10: break
                de=-2*self.spatial*u**-5
                jac=self.L+sp.diags((-7*a2*u**-8/8+2*np.pi*G*(de*u**5+5*e*u**4))[self.idx])
                step=splu(jac.tocsc()).solve(-residual)
                alpha=1.
                while np.min(u[self.idx]+alpha*step)<=0: alpha/=2
                u[self.idx]+=alpha*step
            else: raise RuntimeError('Hamiltonian Newton failed')
            wnew=self.solve_vector(u)
            change=np.max(abs(wnew-w)); w=wnew
            anew=self.longitudinal(w); a2=np.sum(anew*anew,axis=(0,1))
            h=(self.lap@(u-1)+a2*u**-7/8+2*np.pi*G*self.energy(u)*u**5)[self.idx]
            if change<1e-12 and np.max(abs(h))<1e-9: break
        else: raise RuntimeError('coupled iteration failed')
        return u,w,iteration+1

    def diagnostics(self,u,w):
        a=self.longitudinal(w); a2=np.sum(a*a,axis=(0,1))
        # Physical Hamiltonian, recomputed from R3 and KijKij.
        R3=-8*u**-5*(self.lap@u)
        H=R3-u**-12*a2-16*np.pi*G*self.energy(u)
        # Separate divergence-of-tensor route (not expanded vector solve).
        divA=np.array([sum(self.grad[j]@a[i,j] for j in range(3)) for i in range(3)])
        Mi=u**-10*divA-8*np.pi*G*u**-4*self.S
        # Compare on a fixed radial bulk fraction to exclude one-sided edge stencils.
        bulk=(self.r>=RMIN+.25*(RMAX-RMIN))&(self.r<=RMAX-.25*(RMAX-RMIN))
        weights=(self.r*self.q*u**6)[bulk]
        rms=lambda v:float(np.sqrt(np.sum(v[bulk]**2*weights)/np.sum(weights)))
        vector_res=(self.vector@w[:,self.idx].ravel()-(8*np.pi*G*u[self.idx]**6*self.S[:,self.idx]).ravel()).reshape(3,-1)
        return {'resolution':list(self.shape),'all_finite':bool(np.all(np.isfinite(u)) and np.all(np.isfinite(w))),
            'psi_min':float(u.min()),'psi_max':float(u.max()),'H_physical_max_abs_interior':float(np.max(abs(H[self.idx]))),
            'H_physical_bulk_rms':rms(H),'expanded_momentum_solver_max_abs':float(np.max(abs(vector_res))),
            'independent_M_cartesian_bulk_rms':[rms(x) for x in Mi],
            'independent_M_cartesian_max_abs_interior':[float(np.max(abs(x[self.idx]))) for x in Mi],
            'independent_M_vector_bulk_rms':float(np.sqrt(sum(rms(x)**2 for x in Mi))),
            'A_trace_max_abs':float(np.max(abs(np.trace(a,axis1=0,axis2=1)))),
            'boundary_psi_error':float(np.max(abs(u.reshape(self.shape)[[0,-1]]-1))),
            'boundary_W_error':float(np.max(abs(w.reshape((3,)+self.shape)[:,[0,-1]])))}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--levels',type=int,default=3); args=parser.parse_args()
    if args.levels not in (1,2,3): parser.error('--levels must be 1, 2, or 3')
    rows=[]
    for level in range(args.levels):
        nr=8*2**level+1; nt=8*2**level
        print(f'Solving {nr} x {nt} x {nt}',flush=True)
        grid=Grid(nr,nt,nt); u,w,iterations=grid.solve(); row=grid.diagnostics(u,w); row['coupled_iterations']=iterations
        rows.append(row); print(json.dumps(row),flush=True)
        np.savez_compressed(Path(__file__).with_name(f'toroidal_coupled_fields_{nr}_{nt}.npz'),psi=u.reshape(grid.shape),W_cartesian=w.reshape((3,)+grid.shape),resolution=grid.shape)
    orders=[float(np.log2(rows[i]['independent_M_vector_bulk_rms']/rows[i+1]['independent_M_vector_bulk_rms'])) for i in range(len(rows)-1)]
    gate={'finite_positive_solution':all(x['all_finite'] and x['psi_min']>0 for x in rows),
          'physical_hamiltonian_below_1e-7':all(x['H_physical_max_abs_interior']<1e-7 for x in rows),
          'expanded_momentum_below_1e-9':all(x['expanded_momentum_solver_max_abs']<1e-9 for x in rows),
          'boundary_conditions_exact':all(x['boundary_psi_error']==0 and x['boundary_W_error']==0 for x in rows),
          'three_resolution_independent_momentum_convergence':len(rows)==3 and all(p>1.5 for p in orders)}
    out={'checkpoint_name':'Coupled nonlinear toroidal-shell initial-data benchmark','status':'PASS' if all(gate.values()) else 'INCOMPLETE',
      'bulk_radial_interval':[0.10,0.20],
      'raw_field_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('toroidal_coupled_fields_*.npz')},
      'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
      'solver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'frozen_baseline_commit':'2ff3b7ecbd01ad8e2aede2823864d6627275ab43',
      'assumptions':{'G':G,'gamma':'psi^4 times Euclidean toroidal metric','K':0,'A_contravariant':'psi^-10 (L W)','radial_boundary':'psi=1, W=0 at both shell radii','angular_boundary':'periodic','matter':'fixed Run38 fields, covector gradients and normal velocities; E recomputed with physical inverse metric','symmetry':'Run38 stress is axisymmetric; Cartesian W varies with phi'},
      'rows':rows,'independent_momentum_orders':orders,'gate':gate,
      'scope':'finite toroidal shell benchmark; independent momentum residual is truncation error and must converge; no exterior matching, general nonaxisymmetric validation, evolution or physical DEE proof'}
    Path(__file__).with_name('toroidal_helical_coupled_constraints.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    if out['status']!='PASS': raise SystemExit(1)
if __name__=='__main__': main()
