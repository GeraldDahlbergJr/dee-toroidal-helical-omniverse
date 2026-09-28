"""Short-duration reduced coupled Einstein-DEE evolution.

This executable evolves the *physical* 3-metric, extrinsic curvature, and all
three DEE scalar fields together in the local cylindrical-annulus reduction.
It is a genuine coupled ADM Einstein-matter calculation in the stated reduced
symmetry sector; it is NOT yet the full 3-D BSSN/moving-puncture gate described
in docs/COUPLED_EINSTEIN_DEE_EVOLUTION_FORMULATION.md.

Coordinates are (r,phi,z), all geometric coefficients depend only on r, and
Theta/Psi retain fixed linear phi/z gradients while their radial/time profiles
are evolved.  Gauge for this first short cross-check: alpha=1, beta^i=0.

The code starts from the already constraint-solved dynamical DEE checkpoint,
uses the complete action-derived stress tensor, evolves with RK4, and writes a
machine-readable refinement checkpoint.  PASS here means only short-duration
survival of this reduced coupled ADM configuration.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from solve_dynamical_dee_constraints import (
    solve as solve_initial,
    RMIN, RMAX, G, RHO2,
    K_THETA, K_PSI, PI_THETA, PI_PSI,
)
from dee_stress_energy_3p1 import projections, coeffs, potential

M_THETA = 0.0
M_PSI = 0.0
LAMV = 4.0
VEV = 1.0
A = 0.30
B = 0.20
C = 0.10


def deriv1(u, h):
    out = np.empty_like(u)
    out[1:-1] = (u[2:] - u[:-2]) / (2*h)
    out[0] = (u[1] - u[0]) / h
    out[-1] = (u[-1] - u[-2]) / h
    return out


def deriv2(u, h):
    out = np.empty_like(u)
    out[1:-1] = (u[2:] - 2*u[1:-1] + u[:-2]) / h**2
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def spatial_geometry(gamma, h):
    """Return inverse metric, Christoffels, Ricci tensor and scalar.

    gamma has shape (n,3,3) in coordinate basis (r,phi,z), with only radial
    coordinate dependence.  The r^2 cylindrical factor is contained directly
    in gamma_phiphi, so ordinary coordinate Christoffel formulas apply.
    """
    n = len(gamma)
    inv = np.linalg.inv(gamma)
    gp = np.empty_like(gamma)
    for i in range(3):
        for j in range(3):
            gp[:,i,j] = deriv1(gamma[:,i,j], h)

    # Gamma^i_jk
    Chr = np.zeros((n,3,3,3))
    for p in range(n):
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    s = 0.0
                    for l in range(3):
                        dj = gp[p,l,k] if j == 0 else 0.0
                        dk = gp[p,l,j] if k == 0 else 0.0
                        dl = gp[p,j,k] if l == 0 else 0.0
                        s += inv[p,i,l] * (dj + dk - dl)
                    Chr[p,i,j,k] = 0.5*s

    dChr = np.empty_like(Chr)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                dChr[:,i,j,k] = deriv1(Chr[:,i,j,k], h)

    Ric = np.zeros((n,3,3))
    for p in range(n):
        for i in range(3):
            for j in range(3):
                val = dChr[p,0,i,j]
                if j == 0:
                    val -= sum(dChr[p,k,i,k] for k in range(3))
                for k in range(3):
                    tr = sum(Chr[p,l,k,l] for l in range(3))
                    val += Chr[p,k,i,j] * tr
                    for l in range(3):
                        val -= Chr[p,l,i,k] * Chr[p,k,j,l]
                Ric[p,i,j] = val
    R = np.einsum("nij,nij->n", inv, Ric)
    return inv, Chr, Ric, R


def field_metric_and_christoffel(rho):
    kt, kp, la = coeffs(float(rho))
    F = np.array([[1.,0.,0.],[0.,kt,la],[0.,la,kp]])
    Finv = np.linalg.inv(F)
    dF = np.array([[0.,0.,0.],
                   [0.,2*A*rho,2*C*rho],
                   [0.,2*C*rho,2*B*rho]])
    Cfield = np.zeros((3,3,3))
    # only derivative wrt field coordinate 0=rho is nonzero
    for aa in range(3):
        for bb in range(3):
            for cc in range(3):
                s = 0.0
                for dd in range(3):
                    db = dF[dd,cc] if bb == 0 else 0.0
                    dc = dF[dd,bb] if cc == 0 else 0.0
                    ddterm = dF[bb,cc] if dd == 0 else 0.0
                    s += Finv[aa,dd]*(db + dc - ddterm)
                Cfield[aa,bb,cc] = 0.5*s
    return F, Finv, Cfield


def matter_sources(fields, Pi, gamma, inv, h):
    n = len(fields)
    D = np.zeros((n,3,3))
    for a in range(3):
        D[:,a,0] = deriv1(fields[:,a], h)
    D[:,1,1] = M_THETA
    D[:,2,1] = M_PSI
    D[:,1,2] = K_THETA
    D[:,2,2] = K_PSI

    E = np.empty(n); Si = np.empty((n,3)); Sij = np.empty((n,3,3)); S = np.empty(n)
    for p in range(n):
        E[p],Si[p],Sij[p],S[p] = projections(fields[p,0], Pi[p], D[p], inv[p])
    return D,E,Si,Sij,S


def scalar_laplacians(fields, gamma, inv, h):
    n = len(fields)
    det = np.linalg.det(gamma)
    root = np.sqrt(det)
    D = np.zeros((n,3,3))
    for a in range(3):
        D[:,a,0] = deriv1(fields[:,a], h)
    D[:,1,1] = M_THETA; D[:,2,1] = M_PSI
    D[:,1,2] = K_THETA; D[:,2,2] = K_PSI
    laps = np.zeros((n,3))
    for a in range(3):
        flux = np.zeros(n)
        for j in range(3):
            flux += root * inv[:,0,j] * D[:,a,j]
        laps[:,a] = deriv1(flux,h)/root
    return D,laps


def constraints(gamma,Kij,fields,Pi,h):
    inv,Chr,Ric,R = spatial_geometry(gamma,h)
    D,E,Si,Sij,S = matter_sources(fields,Pi,gamma,inv,h)
    Ktr = np.einsum("nij,nij->n",inv,Kij)
    Kup = np.einsum("nik,njl,nkl->nij",inv,inv,Kij)
    K2 = np.einsum("nij,nij->n",Kij,Kup)
    H = R + Ktr*Ktr - K2 - 16*np.pi*G*E

    # M^i = D_j(K^ij-gamma^ij K)-8piG S^i
    Aup = Kup - inv*Ktr[:,None,None]
    dA = np.empty_like(Aup)
    for i in range(3):
        for j in range(3):
            dA[:,i,j] = deriv1(Aup[:,i,j],h)
    M = np.zeros((len(gamma),3))
    Sup = np.einsum("nij,nj->ni",inv,Si)
    for p in range(len(gamma)):
        for i in range(3):
            div = dA[p,i,0]
            for j in range(3):
                for k in range(3):
                    div += Chr[p,i,j,k]*Aup[p,k,j] + Chr[p,j,j,k]*Aup[p,i,k]
            M[p,i] = div - 8*np.pi*G*Sup[p,i]
    return H,M,E,Si,Sij,S,R


def rhs(state,h):
    gamma,Kij,fields,Pi = state
    inv,Chr,Ric,R = spatial_geometry(gamma,h)
    D,E,Si,Sij,S = matter_sources(fields,Pi,gamma,inv,h)
    Ktr = np.einsum("nij,nij->n",inv,Kij)
    Ki_up = np.einsum("nik,nkj->nij",inv,Kij)
    KK = np.einsum("nik,nkj->nij",Kij,Ki_up)

    gdot = -2.0*Kij
    Kdot = Ric + Ktr[:,None,None]*Kij - 2.0*KK
    Kdot += 4*np.pi*G*((S-E)[:,None,None]*gamma - 2.0*Sij)

    D,laps = scalar_laplacians(fields,gamma,inv,h)
    fdot = Pi.copy()
    Pdot = np.zeros_like(Pi)
    # spacetime contraction dphi_B.dphi_C = -Pi_B Pi_C + gamma^ij D_iB D_jC
    for p in range(len(fields)):
        _,Finv,Cfield = field_metric_and_christoffel(fields[p,0])
        contraction = -np.outer(Pi[p],Pi[p]) + np.einsum("bi,ij,cj->bc",D[p],inv[p],D[p])
        gradV = np.array([LAMV*fields[p,0]*(fields[p,0]**2-VEV**2),0.,0.])
        Pdot[p] = Ktr[p]*Pi[p] + laps[p]
        for aa in range(3):
            Pdot[p,aa] += np.einsum("bc,bc->",Cfield[aa],contraction)
        Pdot[p] -= Finv @ gradV
    return gdot,Kdot,fdot,Pdot


def add_state(a,b,scale):
    return tuple(x + scale*y for x,y in zip(a,b))


def rk4_step(state,dt,h,initial_state):
    k1=rhs(state,h)
    k2=rhs(add_state(state,k1,0.5*dt),h)
    k3=rhs(add_state(state,k2,0.5*dt),h)
    k4=rhs(add_state(state,k3,dt),h)
    out=tuple(x + dt*(a+2*b+2*c+d)/6 for x,a,b,c,d in zip(state,k1,k2,k3,k4))
    # Freeze two-point boundary buffer to initial data; diagnostics exclude a wider zone.
    out=[x.copy() for x in out]
    for q in range(4):
        out[q][:2] = initial_state[q][:2]
        out[q][-2:] = initial_state[q][-2:]
    return tuple(out)


def initial_state(n):
    sol=solve_initial()
    r=np.linspace(RMIN,RMAX,n)
    psi,_,W,dW=sol.sol(r)
    gamma=np.zeros((n,3,3))
    gamma[:,0,0]=psi**4
    gamma[:,1,1]=psi**4*r*r
    gamma[:,2,2]=psi**4
    Kij=np.zeros_like(gamma)
    Krz=psi**-2*dW
    Kij[:,0,2]=Krz; Kij[:,2,0]=Krz
    fields=np.zeros((n,3))
    fields[:,0]=np.sqrt(RHO2)
    Pi=np.zeros((n,3))
    Pi[:,1]=PI_THETA; Pi[:,2]=PI_PSI
    return r,(gamma,Kij,fields,Pi)


def norms(H,M,trim=8):
    sl=slice(trim,-trim)
    hm=float(np.max(np.abs(H[sl])))
    hr=float(np.sqrt(np.mean(H[sl]**2)))
    mm=float(np.max(np.abs(M[sl])))
    mr=float(np.sqrt(np.mean(M[sl]**2)))
    return hm,hr,mm,mr


def run(n,t_end,cfl):
    r,state=initial_state(n); initial=tuple(x.copy() for x in state)
    h=r[1]-r[0]; dt=cfl*h; steps=int(np.ceil(t_end/dt)); dt=t_end/steps
    # HARD t=0 interface gate.  The evolved state originates in the reduced
    # cylindrical conformal solver, so timestep 1 is forbidden unless its
    # independently validated constraint evaluator is reproduced first.
    sol=solve_initial()
    from solve_dynamical_dee_constraints import independent_residuals
    ref_h,ref_hr,ref_m,ref_mr=independent_residuals(sol,n)
    # Evaluate the identical reduced conformal operators from the fields that
    # were handed to the evolution code (not by reusing the BVP derivatives).
    psi=np.power(state[0][:,0,0],0.25)
    Krz=state[1][:,0,2]
    dW=psi**2*Krz
    ri=r[1:-1]; p=psi[1:-1]
    dp=(psi[2:]-psi[:-2])/(2*h)
    ddp=(psi[2:]-2*psi[1:-1]+psi[:-2])/h**2
    ddw=(dW[2:]-2*dW[1:-1]+dW[:-2])/h**2
    lap_p=ddp+dp/ri
    lap_w=ddw+dW[1:-1]/ri
    abar2=2*dW[1:-1]**2
    from solve_dynamical_dee_constraints import energy, S_Z
    Hred=-8*p**-5*lap_p-p**-12*abar2-16*np.pi*G*energy(p)
    Mred=p**-10*lap_w-8*np.pi*G*p**-4*S_Z
    evo_h=float(np.max(np.abs(Hred))); evo_m=float(np.max(np.abs(Mred)))
    # Same second-order stencil and same physical conventions should agree
    # closely with the independent evaluator; tolerance is deliberately tight.
    rel_h=abs(evo_h-ref_h)/max(ref_h,1e-30)
    rel_m=abs(evo_m-ref_m)/max(ref_m,1e-30)
    interface_pass=(rel_h < 0.03 and rel_m < 0.03)
    if not interface_pass:
        return {"points":n,"interface_pass":False,
                "reference_hamiltonian_max_abs":ref_h,
                "evolution_t0_hamiltonian_max_abs":evo_h,
                "reference_momentum_max_abs":ref_m,
                "evolution_t0_momentum_max_abs":evo_m,
                "hamiltonian_relative_mismatch":rel_h,
                "momentum_relative_mismatch":rel_m,
                "finite":True,"aborted_before_timestep_1":True}
    # Keep the generic ADM evaluator as a diagnostic, but do not confuse it
    # with the validated reduced-conformal t=0 gate.
    H0,M0,*_=constraints(*state,h); n0=norms(H0,M0)
    maxH=n0[0]; maxM=n0[2]; rho_min=float(np.min(state[2][:,0]))
    min_det=float(np.min(np.linalg.det(state[0])))
    for _ in range(steps):
        state=rk4_step(state,dt,h,initial)
        H,M,*_=constraints(*state,h)
        hm,hr,mm,mr=norms(H,M)
        maxH=max(maxH,hm); maxM=max(maxM,mm)
        rho_min=min(rho_min,float(np.min(state[2][8:-8,0])))
        min_det=min(min_det,float(np.min(np.linalg.det(state[0][8:-8]))))
        if not np.isfinite(maxH+maxM+rho_min+min_det):
            break
    H,M,*_=constraints(*state,h); nf=norms(H,M)
    return {
        "points":n,"dt":dt,"steps":steps,"interface_pass":True,"aborted_before_timestep_1":False,
        "reference_hamiltonian_max_abs":ref_h,"evolution_t0_hamiltonian_max_abs":evo_h,
        "reference_momentum_max_abs":ref_m,"evolution_t0_momentum_max_abs":evo_m,
        "hamiltonian_relative_mismatch":rel_h,"momentum_relative_mismatch":rel_m,
        "initial_hamiltonian_max_abs":n0[0],"initial_momentum_max_abs":n0[2],
        "final_hamiltonian_max_abs":nf[0],"final_hamiltonian_rms":nf[1],
        "final_momentum_max_abs":nf[2],"final_momentum_rms":nf[3],
        "max_hamiltonian_over_run":maxH,"max_momentum_over_run":maxM,
        "rho_min":rho_min,"min_det_gamma":min_det,
        "finite":bool(np.isfinite(maxH) and np.isfinite(maxM) and np.isfinite(rho_min) and np.isfinite(min_det))
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--resolutions",nargs="+",type=int,default=[201,401,801])
    ap.add_argument("--t-end",type=float,default=.25)
    ap.add_argument("--cfl",type=float,default=.10)
    ap.add_argument("--output",default="coupled_einstein_dee_checkpoint.json")
    args=ap.parse_args()
    rows=[run(n,args.t_end,args.cfl) for n in args.resolutions]
    # Convergence of the t=0 interface itself is always recorded.
    for key in ("reference_hamiltonian","reference_momentum","evolution_t0_hamiltonian","evolution_t0_momentum"):
        prev=None
        for row in rows:
            cur=row[f"{key}_max_abs"]
            row[f"{key}_observed_order"]=None if prev is None else float(np.log(prev/cur)/np.log(2))
            prev=cur
    if not all(x.get("interface_pass",False) for x in rows):
        out={"checkpoint_name":"Coupled Einstein-DEE t=0 interface gate",
             "status":"INCONCLUSIVE","scope":"pre-evolution cylindrical conformal consistency gate",
             "runs":rows,
             "gate":{"t0_interface_consistent":False},
             "next_gate":"No timestep permitted until this passes."}
        Path(args.output).write_text(json.dumps(out,indent=2)+"\\n")
        print(json.dumps(out,indent=2))
        raise SystemExit(3)
    # Only after the hard interface gate passes do we assess evolved residuals.
    for key in ("hamiltonian","momentum"):
        prev=None
        for row in rows:
            cur=row[f"final_{key}_max_abs"]
            row[f"{key}_observed_order"]=None if prev is None else float(np.log(prev/cur)/np.log(2))
            prev=cur
    gate={
        "t0_interface_consistent":all(x.get("interface_pass",False) for x in rows),
        "all_runs_finite":all(x["finite"] for x in rows),
        "metric_positive":all(x["min_det_gamma"]>0 for x in rows),
        "finite_branch_positive":all(x["rho_min"]>0 for x in rows),
        "hamiltonian_refines":all(rows[i]["final_hamiltonian_max_abs"]<rows[i-1]["final_hamiltonian_max_abs"] for i in range(1,len(rows))),
        "momentum_refines":all(rows[i]["final_momentum_max_abs"]<rows[i-1]["final_momentum_max_abs"] for i in range(1,len(rows))),
    }
    if len(rows)>=3:
        gate["hamiltonian_final_order_positive"]=rows[-1]["hamiltonian_observed_order"]>0
        gate["momentum_final_order_positive"]=rows[-1]["momentum_observed_order"]>0
    status="PASS" if all(gate.values()) else ("FAIL" if all(x["finite"] for x in rows) else "INCONCLUSIVE")
    out={
        "checkpoint_name":"Short-duration reduced coupled Einstein-DEE ADM evolution",
        "status":status,
        "scope":"1D cylindrical-annulus ADM evolution with alpha=1 beta=0; matter and geometry both evolved",
        "not_claimed":"full 3-D BSSN/moving-puncture evolution, global stability, or experimental validation",
        "t_end":args.t_end,"cfl":args.cfl,"G":G,
        "runs":rows,"gate":gate
    }
    Path(args.output).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
    # A numerical FAIL is a scientific result, not an infrastructure failure.
    # Exit nonzero only for non-finite/inconclusive execution.
    if status=="INCONCLUSIVE":
        raise SystemExit(3)

if __name__=="__main__":
    main()
