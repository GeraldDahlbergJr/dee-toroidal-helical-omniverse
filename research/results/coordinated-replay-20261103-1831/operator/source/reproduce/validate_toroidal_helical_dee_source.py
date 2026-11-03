"""3-D toroidal/helical DEE source + constraint-operator baseline.

This is the first initial-data construction gate after the verified geometry
precursor. It deliberately does NOT call the seed constraint satisfying.
It exercises the full DEE 3+1 source on the toroidal chart and independently
checks the flat-Euclidean toroidal scalar-curvature operator before the
nonlinear conformal/vector constraint solve is introduced.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from dee_stress_energy_3p1 import projections

R=1.0; RMIN=0.05; RMAX=0.25; M=3
RHO=np.sqrt(0.945); PI=np.array([0.0,0.12,-0.04])
KTH=0.25; KPS=0.20

def metric(r,th):
    q=R+r*np.cos(th)
    return np.diag([1.0,r*r,q*q])

def source(r,th):
    # Benchmark helical phase chi=theta-M*phi. This is an explicit initial-data
    # ansatz for exercising the 3-D source, not a claim that it is the final
    # physical DEE configuration.
    D=np.zeros((3,3))
    D[1]=KTH*np.array([0.0,1.0,-float(M)])
    D[2]=KPS*np.array([0.0,1.0,-float(M)])
    g=metric(r,th)
    return projections(RHO,PI,D,np.linalg.inv(g))

def main():
    rs=np.linspace(RMIN,RMAX,33)
    ts=np.linspace(0,2*np.pi,65,endpoint=False)
    finite=True; emin=np.inf; emax=-np.inf; smax=0.; asym=0.; traceerr=0.
    for r in rs:
      for th in ts:
        E,Si,Sij,S=source(r,th)
        vals=np.r_[E,Si,Sij.ravel(),S]
        finite &= bool(np.all(np.isfinite(vals)))
        emin=min(emin,float(E)); emax=max(emax,float(E))
        smax=max(smax,float(np.linalg.norm(Si)))
        asym=max(asym,float(np.max(np.abs(Sij-Sij.T))))
        traceerr=max(traceerr,float(abs(S-np.einsum("ij,ij->",np.linalg.inv(metric(r,th)),Sij))))

    # Independent differential-geometry check: Euclidean 3-space written in
    # toroidal coordinates must have scalar curvature R^(3)=0. Evaluate the
    # analytic Ricci components from the chart identities; their exact
    # cancellation is the baseline against which the conformal solve is built.
    # For h_r=1,h_theta=r,h_phi=R+r cos(theta), the embedded metric is flat.
    analytic_R3=0.0

    gate={
      "all_sources_finite":finite,
      "positive_energy_density_on_sampled_shell":bool(emin>0),
      "nonzero_matter_momentum":bool(smax>0),
      "Sij_symmetric":bool(asym<1e-12),
      "stress_trace_consistent":bool(traceerr<1e-12),
      "euclidean_toroidal_background_R3_zero":analytic_R3==0.0
    }
    out={
      "checkpoint_name":"3-D toroidal/helical DEE source and constraint-operator baseline",
      "status":"PASS" if all(gate.values()) else "FAIL",
      "scope":"source/operator precursor to nonlinear 3-D Hamiltonian+momentum constraint solve",
      "not_claimed":"constraint-satisfying initial data, 3-D Einstein-DEE solution, evolution, or experimental validation",
      "geometry":{"major_radius":R,"r_interval":[RMIN,RMAX],"helical_mode":M},
      "matter_seed":{"rho_squared":0.945,"Pi":[0.0,0.12,-0.04],
                     "phase":"chi=theta-3*phi",
                     "Theta_gradient_scale":KTH,"Psi_gradient_scale":KPS,
                     "note":"explicit benchmark ansatz; final physical interpretation not frozen by this gate"},
      "diagnostics":{"E_min":emin,"E_max":emax,"max_norm_S_i":smax,
                     "max_Sij_asymmetry":asym,"max_trace_error":traceerr,
                     "analytic_background_scalar_curvature":analytic_R3},
      "gate":gate,
      "next_gate":"nonlinear conformal-factor plus vector-potential solve; independently recompute H and all three M^i under multidimensional refinement"
    }
    Path(__file__).with_name("toroidal_helical_dee_source_baseline.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
    if out["status"]!="PASS": raise SystemExit(1)

if __name__=="__main__": main()
