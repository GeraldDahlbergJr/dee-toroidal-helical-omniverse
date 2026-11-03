"""Coupled conformal Hamiltonian/vector benchmark on the Run 38 shell.

Cartesian components of W on a toroidal grid avoid component-connection
ambiguity. gamma=psi^4 delta, K=0, A^ij=psi^-10 (LW)^ij. Fixed physical
normal field velocities and coordinate gradients; E is recomputed with gamma.
The frozen seed has no phi-dependent stress despite its helical phase.
This is a 3-coordinate operator benchmark, not generic 3-D matter validation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import newton_krylov, NoConvergence
from scipy.sparse.linalg import LinearOperator
from scipy.fft import dst, idst, fftn, ifftn
from validate_toroidal_helical_dee_source import R,RMIN,RMAX,source,RHO,PI
from dee_stress_energy_3p1 import coeffs,potential
G=1e-3

def solve(nr,nt,nph):
    shape=(nr,nt,nph); interior=(nr-2,nt,nph)
    r=np.linspace(RMIN,RMAX,nr)[:,None,None]
    t=np.arange(nt)[None,:,None]*2*np.pi/nt
    p=np.arange(nph)[None,None,:]*2*np.pi/nph
    dr=(RMAX-RMIN)/(nr-1); dt=2*np.pi/nt; dp=2*np.pi/nph
    q=R+r*np.cos(t)
    er=np.array(np.broadcast_arrays(np.cos(t)*np.cos(p),np.cos(t)*np.sin(p),np.sin(t)+0*p+0*r))
    et=np.array(np.broadcast_arrays(-np.sin(t)*np.cos(p),-np.sin(t)*np.sin(p),np.cos(t)+0*p+0*r))
    ep=np.array(np.broadcast_arrays(-np.sin(p)+0*r+0*t,np.cos(p)+0*r+0*t,0*r+0*t+0*p))
    E0=np.empty(shape); Scart=np.empty((3,)+shape)
    for i in range(nr):
      for j in range(nt):
        e,s,_,_=source(float(r[i,0,0]),float(t[0,j,0]))
        E0[i,j]=e
        Scart[:,i,j]=er[:,i,j]*s[0]+et[:,i,j]*s[1]/r[i,0,0]+ep[:,i,j]*s[2]/q[i,j,0]
    kt,kp,la=coeffs(RHO); kin=np.array([[1,0,0],[0,kt,la],[0,la,kp]])
    ec=.5*float(PI@kin@PI)+potential(RHO); eg=E0-ec
    def deriv(u):
      return np.gradient(u,dr,axis=0,edge_order=2),(np.roll(u,-1,1)-np.roll(u,1,1))/(2*dt),(np.roll(u,-1,2)-np.roll(u,1,2))/(2*dp)
    def grad(u):
      ur,ut,up=deriv(u)
      return er*ur+et*(ut/r)+ep*(up/q)
    def lap(u):
      out=np.zeros(shape)
      out[1:-1]=(u[2:]-2*u[1:-1]+u[:-2])/dr**2+(1/r[1:-1]+np.cos(t)/q[1:-1])*(u[2:]-u[:-2])/(2*dr)
      out+=(np.roll(u,-1,1)-2*u+np.roll(u,1,1))/(dt*dt*r*r)-np.sin(t)/(r*q)*(np.roll(u,-1,1)-np.roll(u,1,1))/(2*dt)
      out+=(np.roll(u,-1,2)-2*u+np.roll(u,1,2))/(dp*dp*q*q)
      return out
    # Manufactured Cartesian polynomial fields exercise all angular directions.
    xyz=np.array(np.broadcast_arrays(q*np.cos(p),q*np.sin(p),r*np.sin(t)+0*p))
    trial=xyz**2
    trial_div=sum(grad(trial[i])[i] for i in range(3))
    trial_vec=np.array([lap(v) for v in trial])+grad(trial_div)/3
    band=(r[:,0,0]>=.10-1e-12)&(r[:,0,0]<=.20+1e-12)
    poly_error=float(np.sqrt(np.mean((trial_vec[:,band]-8/3)**2)))
    lap_error=float(np.sqrt(np.mean((lap(np.sum(xyz**2,axis=0))[band]-6)**2)))

    def unpack(x):
      z=x.reshape((4,)+interior); u=np.ones(shape); u[1:-1]=z[0]
      w=np.zeros((3,)+shape); w[:,1:-1]=z[1:]
      return u,w
    def tensors(w):
      dw=np.array([grad(v) for v in w]); div=np.einsum('ii...->...',dw)
      aa=dw+np.swapaxes(dw,0,1)-2/3*np.eye(3)[:,:,None,None,None]*div
      return div,aa
    def residual(x):
      u,w=unpack(x); div,aa=tensors(w)
      en=ec+eg*u**-4; aa2=np.einsum('ij...,ij...->...',aa,aa)
      h=lap(u)+aa2*u**-7/8+2*np.pi*G*en*u**5
      mom=np.array([lap(v) for v in w])+grad(div)/3-8*np.pi*G*u**6*Scart
      return np.concatenate((h[None,1:-1],mom[:,1:-1]),axis=0).ravel()
    # Constant-coefficient inverse Laplacian: sine radial/FFT angular.
    lr=-4*np.sin(np.pi*np.arange(1,nr-1)/(2*(nr-1)))**2/dr**2
    lt=-4*np.sin(np.pi*np.arange(nt)/nt)**2/(dt*dt*((RMIN+RMAX)/2)**2)
    lp=-4*np.sin(np.pi*np.arange(nph)/nph)**2/(dp*dp*R*R)
    eig=lr[:,None,None]+lt[None,:,None]+lp[None,None,:]
    def inverse(v):
      a=v.reshape((4,)+interior)
      tr=fftn(dst(a,type=1,axis=1,norm='ortho'),axes=(2,3))
      return idst(ifftn(tr/eig[None],axes=(2,3)).real,type=1,axis=1,norm='ortho').ravel()
    size=4*np.prod(interior); pre=LinearOperator((size,size),matvec=inverse)
    x=np.zeros((4,)+interior); x[0]=1
    try:
      x=newton_krylov(residual,x.ravel(),f_tol=1e-9,maxiter=80,inner_M=pre)
      converged=True
    except NoConvergence as e:
      x=e.args[0]; converged=False
    u,w=unpack(x); div,aa=tensors(w); f=residual(x).reshape((4,)+interior)
    h=-8*u[1:-1]**-5*f[0]
    # Independently differentiate the constructed LW tensor; unlike the solve,
    # this route does not use lap(W)+grad(div W)/3. Report bulk and full shell.
    diva=np.zeros((3,)+shape)
    for i in range(3):
      for j in range(3): diva[i]+=grad(aa[i,j])[j]
    mind=diva-8*np.pi*G*u**6*Scart
    sl=(slice(None),slice(2,-2),slice(None),slice(None))
    # Independent Hamiltonian: conservative face flux.
    rf=(r[1:]+r[:-1])/2; qf=R+rf*np.cos(t)
    fr=rf*qf*(u[1:]-u[:-1])/dr
    lf=(fr[1:]-fr[:-1])/dr/(r[1:-1]*q[1:-1])
    qt=R+r*np.cos(t+dt/2); ft=qt/r*(np.roll(u,-1,1)-u)/dt
    lf+=((ft-np.roll(ft,1,1))/dt/(r*q))[1:-1]
    lf+=((np.roll(u,-1,2)-2*u+np.roll(u,1,2))/(dp*dp*q*q))[1:-1]
    aa2=np.einsum('ij...,ij...->...',aa,aa)
    hind=-8*u[1:-1]**-5*lf-aa2[1:-1]*u[1:-1]**-12-16*np.pi*G*(ec+eg[1:-1]*u[1:-1]**-4)
    return {'resolution':[nr,nt,nph],'solver_converged':converged,'all_finite':bool(np.all(np.isfinite(x))),
      'psi_min':float(u.min()),'psi_max':float(u.max()),'W_max_abs':float(np.max(abs(w))),
      'physical_ADM_H_max_abs':float(np.max(abs(h))),
      'conformal_momentum_solve_max_abs':[float(np.max(abs(f[i]))) for i in range(1,4)],
      'independent_flux_H_rms':float(np.sqrt(np.mean(hind**2))),
      'manufactured_vector_operator_rms':poly_error,
      'manufactured_scalar_operator_rms':lap_error,
      'independent_momentum_matched_band_rms':[float(np.sqrt(np.mean((u**-10*mind[i])[band]**2))) for i in range(3)],
      'independent_tensor_momentum_bulk_rms':[float(np.sqrt(np.mean(mind[sl][i]**2))) for i in range(3)],
      'independent_tensor_momentum_full_rms':[float(np.sqrt(np.mean(mind[i,1:-1]**2))) for i in range(3)]}

def main():
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--levels',type=int,default=3); args=ap.parse_args()
    rows=[]
    for res in [(9,16,16),(17,32,32),(33,64,64)][:args.levels]:
      row=solve(*res); rows.append(row); print(json.dumps(row),flush=True)
    gate={'all_solver_constraints_small':all(r['solver_converged'] and r['physical_ADM_H_max_abs']<1e-7 and max(r['conformal_momentum_solve_max_abs'])<1e-8 for r in rows),
      'positive_finite_metric':all(r['all_finite'] and r['psi_min']>0 for r in rows),
      'manufactured_operators_refine':len(rows)==3 and all(rows[i+1][key]<rows[i][key] for key in ['manufactured_vector_operator_rms','manufactured_scalar_operator_rms'] for i in range(2)),
      'independent_residuals_refine':len(rows)==3 and all(rows[i+1]['independent_flux_H_rms']<rows[i]['independent_flux_H_rms'] and all(rows[i+1]['independent_momentum_matched_band_rms'][k]<rows[i]['independent_momentum_matched_band_rms'][k] for k in range(3)) for i in range(2))}
    orders={key:[float(np.log2(rows[i][key]/rows[i+1][key])) for i in range(len(rows)-1)] for key in ['manufactured_vector_operator_rms','manufactured_scalar_operator_rms','independent_flux_H_rms']}
    out={'observed_refinement_orders':orders,'checkpoint_name':'Coupled nonlinear toroidal Hamiltonian and three momentum benchmark',
      'status':'SOLVER_PASS_VALIDATION_PENDING' if all(gate.values()) else 'FAIL',
      'source_baseline_commit':'2ff3b7ecbd01ad8e2aede2823864d6627275ab43',
      'assumptions':{'G':G,'K':0,'metric':'psi^4 times Euclidean toroidal metric','A':'psi^-10 LW, Cartesian components','radial_boundary':'psi=1 and W=0 at both shell radii','angular_boundaries':'periodic','matter':'Run 38 fixed Pi and gradients, energy recomputed in physical metric'},
      'rows':rows,'gate':gate,
      'limitations':['Frozen seed stress is phi-independent; generic phi-dependent matter remains untested.','Independent tensor-divergence residual and boundary convergence must be assessed before claiming constraint-satisfying continuum data.','Finite-shell benchmark boundaries and G=0.001 are numerical assumptions.']}
    Path(__file__).with_name('toroidal_helical_full_constraints.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)
    if not all(gate.values()): raise SystemExit(1)
if __name__=='__main__': main()
