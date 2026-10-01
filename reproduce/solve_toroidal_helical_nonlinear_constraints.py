"""Nonlinear 3-D toroidal/helical Hamiltonian-constraint construction gate.

Uses the frozen Run #38 benchmark matter profile on the verified Euclidean
toroidal chart. This first nonlinear gate solves the time-symmetric conformal
Hamiltonian equation with K_ij=0. Because Run #38 has nonzero S_i, this is NOT
yet a full constraint-satisfying Einstein-DEE data set: the momentum residuals
are measured explicitly and are expected to remain nonzero. The purpose is to
validate the nonlinear 3-D Hamiltonian solve before introducing the vector
potential required by the momentum constraints.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import newton_krylov
from validate_toroidal_helical_dee_source import source,R,RMIN,RMAX,M

G=1.0e-3
RESOLUTIONS=((9,16,16),(17,32,32),(33,64,64))

def solve(nr,nt,np_):
    r=np.linspace(RMIN,RMAX,nr); th=np.linspace(0,2*np.pi,nt,endpoint=False)
    ph=np.linspace(0,2*np.pi,np_,endpoint=False)
    dr=r[1]-r[0]; dt=2*np.pi/nt; dp=2*np.pi/np_
    rr=r[:,None,None]; tt=th[None,:,None]
    E=np.empty((nr,nt,np_)); S=np.empty((3,nr,nt,np_))
    for i,x in enumerate(r):
      for j,y in enumerate(th):
        e,s,_,_=source(x,y); E[i,j,:]=e; S[:,i,j,:]=s[:,None]

    # Hold normal field velocities and coordinate gradients fixed, and recompute
    # physical energy with gamma^ij=psi^-4 * background_gamma^ij.
    from dee_stress_energy_3p1 import coeffs, potential
    from validate_toroidal_helical_dee_source import RHO, PI
    kt,kp,la=coeffs(RHO)
    kinetic=np.array([[1.,0.,0.],[0.,kt,la],[0.,la,kp]])
    e_constant=.5*float(PI@kinetic@PI)+potential(RHO)
    e_gradient=E-e_constant
    def physical_energy(u):
      return e_constant+e_gradient*u**-4

    # Dirichlet psi=1 at shell radii; theta/phi periodic. This is an explicit
    # numerical benchmark BC, not yet a frozen physical boundary condition.
    def lap(u):
      q=R+rr*np.cos(tt)
      # scalar Laplacian in orthogonal toroidal coordinates
      radial=np.zeros_like(u)
      radial[1:-1]=(u[2:]-2*u[1:-1]+u[:-2])/dr**2
      radial[1:-1]+=(1/rr[1:-1]+np.cos(tt)/q[1:-1])*(u[2:]-u[:-2])/(2*dr)
      theta=(np.roll(u,-1,1)-2*u+np.roll(u,1,1))/(dt*dt*rr*rr)
      theta-=np.sin(tt)/(rr*q)*(np.roll(u,-1,1)-np.roll(u,1,1))/(2*dt)
      phi=(np.roll(u,-1,2)-2*u+np.roll(u,1,2))/(dp*dp*q*q)
      return radial+theta+phi
    shape=(nr-2,nt,np_)
    def residual(x):
      u=np.ones((nr,nt,np_)); u[1:-1]=x.reshape(shape)
      L=lap(u)
      H=L[1:-1]+2*np.pi*G*physical_energy(u)[1:-1]*u[1:-1]**5
      return H.ravel()
    x0=np.ones(shape).ravel()
    x=newton_krylov(residual,x0,f_tol=2e-9,maxiter=50)
    u=np.ones((nr,nt,np_)); u[1:-1]=np.asarray(x).reshape(shape)
    h=residual(x)
    # Physical ADM Hamiltonian for K_ij=0: R-16*pi*G*E.
    physical_h=-8*u[1:-1]**-5*lap(u)[1:-1]-16*np.pi*G*physical_energy(u)[1:-1]
    # A separate conservative face-flux discretization checks consistency as
    # the mesh is refined; it is not used by the nonlinear solve.
    q=R+rr*np.cos(tt)
    rf=.5*(rr[1:]+rr[:-1]); qf=R+rf*np.cos(tt)
    fr=rf*qf*(u[1:]-u[:-1])/dr
    lflux=(fr[1:]-fr[:-1])/dr/(rr[1:-1]*q[1:-1])
    tf=tt+dt/2; qt=R+rr*np.cos(tf)
    ft=qt/rr*(np.roll(u,-1,1)-u)/dt
    lt=(ft-np.roll(ft,1,1))/dt/(rr*q)
    lp=(np.roll(u,-1,2)-2*u+np.roll(u,1,2))/(dp*dp*q*q)
    independent_h=-8*u[1:-1]**-5*(lflux+lt[1:-1]+lp[1:-1])-16*np.pi*G*physical_energy(u)[1:-1]
    # With A_ij=K=0, M_i=-8 pi G S_i. Record all three components rather
    # than pretending this Hamiltonian-only precursor solves momentum.
    mom=-8*np.pi*G*S[:,1:-1]
    return {"resolution":[nr,nt,np_],"all_finite":bool(np.all(np.isfinite(u))),
      "psi_min":float(u.min()),"psi_max":float(u.max()),
      "hamiltonian_residual_max_abs":float(np.max(np.abs(h))),
      "physical_ADM_H_max_abs":float(np.max(np.abs(physical_h))),
      "independent_flux_H_rms":float(np.sqrt(np.mean(independent_h**2))),
      "hamiltonian_residual_rms":float(np.sqrt(np.mean(h*h))),
      "momentum_residual_max_abs":[float(np.max(np.abs(mom[k]))) for k in range(3)]}

def main():
    rows=[solve(*n) for n in RESOLUTIONS]
    row=rows[-1]
    gate={"finite_solution":row["all_finite"],"positive_conformal_factor":row["psi_min"]>0,
          "hamiltonian_solved":all(r["physical_ADM_H_max_abs"]<1e-7 for r in rows),
          "independent_residual_decreases":all(rows[i+1]["independent_flux_H_rms"]<rows[i]["independent_flux_H_rms"] for i in range(len(rows)-1)),
          "momentum_not_silently_claimed_solved":max(row["momentum_residual_max_abs"])>1e-8}
    out={"checkpoint_name":"Nonlinear 3-D toroidal/helical Hamiltonian initial-data precursor",
      "status":"PASS" if all(gate.values()) else "FAIL",
      "scope":"nonlinear Hamiltonian solve on frozen Run #38 benchmark; momentum constraints explicitly unsolved",
      "assumptions":{"conformal_background":"verified Euclidean toroidal chart","mean_curvature_K":0.0,
        "tracefree_extrinsic_curvature":"zero in this precursor","radial_boundary":"psi=1 Dirichlet benchmark",
        "angular_boundaries":"periodic","matter_energy":"recomputed using inverse conformal metric; Pi and coordinate gradients fixed", "G":G},
      "rows":rows,"gate":gate,
      "not_claimed":"full 3-D constraint-satisfying initial data, physical boundary condition, evolution, or experimental validation",
      "next_gate":"introduce vector potential and solve the three momentum constraints together with Hamiltonian constraint; then independent multidimensional refinement"}
    Path(__file__).with_name("toroidal_helical_nonlinear_hamiltonian.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
    if out["status"]!="PASS": raise SystemExit(1)
if __name__=="__main__": main()
