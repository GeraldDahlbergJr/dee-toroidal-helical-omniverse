"""3-D toroidal/helical geometry precursor gate.

This gate is deliberately kinematic. It validates the Euclidean toroidal chart
and the existing dee_geometry.py helical centerline/path convention before any
3-D Einstein-DEE constraint solve is attempted. It does not alter or claim a
spacetime solution, matter initial data, boundary condition, or evolution law.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

R = 1.0
A = 0.25
M = 3
SAMPLES = 4097

def embed(r, theta, phi):
    q = R + r*np.cos(theta)
    return np.stack((q*np.cos(phi), q*np.sin(phi), r*np.sin(theta)), axis=-1)

def metric_diag(r, theta):
    return np.stack((np.ones_like(theta), r*r*np.ones_like(theta),
                     (R+r*np.cos(theta))**2), axis=-1)

def main():
    # Existing dee_geometry path is r=A, theta=M*phi.
    phi=np.linspace(0.0,2.0*np.pi,SAMPLES)
    theta=M*phi
    xyz=embed(A,theta,phi)
    q=R+A*np.cos(theta)
    analytic_speed2=q*q+(A*M)**2

    # Analytic Cartesian tangent to avoid conflating the gate with FD error.
    dq=-A*M*np.sin(theta)
    dx=dq*np.cos(phi)-q*np.sin(phi)
    dy=dq*np.sin(phi)+q*np.cos(phi)
    dz=A*M*np.cos(theta)
    cart_speed2=dx*dx+dy*dy+dz*dz
    pullback_rel=float(np.max(np.abs(cart_speed2-analytic_speed2)/
                              np.maximum(analytic_speed2,1e-30)))

    # Independently reconstruct the existing public trajectory formula.
    radius=R+A*np.cos(M*phi)
    legacy=np.stack((radius*np.cos(phi),radius*np.sin(phi),
                     A*np.sin(M*phi)),axis=-1)
    path_max=float(np.max(np.abs(xyz-legacy)))

    # Sample the nonsingular toroidal shell, excluding r=0 coordinate axis.
    rs=np.linspace(0.05,A,33)
    th=np.linspace(0.0,2.0*np.pi,257)
    min_det=np.inf; min_eig=np.inf; max_sym=0.0
    for r in rs:
        diag=metric_diag(r,th)
        det=diag[:,0]*diag[:,1]*diag[:,2]
        min_det=min(min_det,float(np.min(det)))
        min_eig=min(min_eig,float(np.min(diag)))
    # Closure and regularity of the embedded helical path.
    closure=float(np.linalg.norm(xyz[-1]-xyz[0]))
    finite=bool(np.all(np.isfinite(xyz)) and np.all(np.isfinite(analytic_speed2)))
    positive=bool(min_det>0.0 and min_eig>0.0 and R>A)
    gate=bool(finite and positive and path_max<1e-13 and
              pullback_rel<1e-13 and closure<1e-12)

    out={
      "scope":"3-D toroidal/helical Euclidean geometry precursor gate",
      "status":"PASS" if gate else "FAIL",
      "not_claimed":"3-D Einstein-DEE constraint solution, spacetime solution, evolution, or experimental validation",
      "parameters":{"major_radius":R,"minor_radius":A,"helical_mode":M,"samples":SAMPLES},
      "chart":{"coordinates":["r","theta","phi"],
               "embedding":["(R+r*cos(theta))*cos(phi)","(R+r*cos(theta))*sin(phi)","r*sin(theta)"],
               "metric_diagonal":["1","r^2","(R+r*cos(theta))^2"],
               "sampled_r_interval":[0.05,A]},
      "helical_path":{"definition":"r=A; theta=M*phi",
                      "max_abs_difference_from_dee_geometry_formula":path_max,
                      "closure_error":closure,
                      "max_relative_pullback_line_element_error":pullback_rel},
      "regularity":{"all_finite":finite,"R_greater_than_A":bool(R>A),
                    "minimum_sampled_metric_determinant":min_det,
                    "minimum_sampled_metric_eigenvalue":min_eig},
      "gate":{"metric_positive_definite_on_sampled_shell":positive,
              "matches_existing_helical_path":bool(path_max<1e-13),
              "pullback_matches_cartesian_tangent":bool(pullback_rel<1e-13),
              "closed_path":bool(closure<1e-12)}
    }
    p=Path(__file__).with_name("toroidal_helical_geometry_gate.json")
    p.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
    if not gate: raise SystemExit(1)

if __name__=="__main__":
    main()
