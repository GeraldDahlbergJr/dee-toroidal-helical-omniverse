#!/usr/bin/env python3
"""Off-line Hamiltonian source perturbation diagnostic; NEVER modifies Run43.

Fixed-geometry first variation only:
H = R3 + K^2 - Kij*Kij - 16*pi*G*rho
delta H = -16*pi*G*delta rho.
The input rho is the Eulerian energy density (not necessarily rest-mass density).
No physical evolution, geometry backreaction, or constraint solve is performed.
"""
import argparse
import json
import math


def source_perturbation(r, phi, z, t, *, amplitude, radial_center,
                        radial_width, m, kz, omega, phase):
    if radial_width <= 0:
        raise ValueError("radial_width must be positive")
    envelope = math.exp(-0.5*((r-radial_center)/radial_width)**2)
    return amplitude*envelope*math.cos(m*phi+kz*z+omega*t+phase)


def evaluate(*, r, phi, z, t, amplitude, radial_center, radial_width,
             m, kz, omega, phase, gravitational_constant):
    if gravitational_constant <= 0:
        raise ValueError("G must be positive")
    vals=(r,phi,z,t,amplitude,radial_center,radial_width,kz,omega,phase,gravitational_constant)
    if not all(math.isfinite(x) for x in vals):
        raise ValueError("all numerical inputs must be finite")
    delta_rho=source_perturbation(r,phi,z,t,amplitude=amplitude,
            radial_center=radial_center,radial_width=radial_width,
            m=m,kz=kz,omega=omega,phase=phase)
    return {"delta_eulerian_energy_density":delta_rho,
            "fixed_geometry_delta_hamiltonian":-16*math.pi*gravitational_constant*delta_rho,
            "helical_phase":"m*phi+kz*z+omega*t+phase",
            "limitations":"Fixed-geometry source-only sensitivity. Not constraint-satisfying initial data. No backreaction or evolution."}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name,default in (("r",1.0),("phi",0.0),("z",0.0),("t",0.0),
                         ("amplitude",1e-6),("radial-center",1.0),
                         ("radial-width",0.2),("kz",1.0),("omega",1.0),
                         ("phase",0.0),("gravitational-constant",1.0)):
        p.add_argument("--"+name,type=float,default=default)
    p.add_argument("--m",type=int,default=2)
    a=vars(p.parse_args())
    a={k.replace("-","_"):v for k,v in a.items()}
    print(json.dumps(evaluate(**a),indent=2,allow_nan=False))


if __name__=="__main__":
    main()
