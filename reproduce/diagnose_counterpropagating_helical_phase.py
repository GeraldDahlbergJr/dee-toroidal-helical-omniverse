#!/usr/bin/env python3
"""Diagnostic-only counter-propagating helical phase experiment.

This is a kinematic toy model, NOT a solution of Einstein/DEE field equations.
No production parameters, matter sources, or validation thresholds are changed.
Uses identical helical spatial phase with opposite temporal phase propagation.
"""
import argparse
import json
import math


def sample(m=2, kz=1.0, omega=1.0, forward_amplitude=1.0,
           reverse_amplitude=1.0, relative_phase=0.0, n=256, t=0.0):
    if n < 4 or not all(math.isfinite(x) for x in
            (kz, omega, forward_amplitude, reverse_amplitude, relative_phase, t)):
        raise ValueError("invalid input")
    # Coordinate s is the spatial helical phase m*phi+kz*z.
    # Opposite signs of omega represent opposite phase motion in s.
    rows = []
    for j in range(n):
        s = 2*math.pi*j/n
        forward = forward_amplitude*math.cos(s-omega*t)
        reverse = reverse_amplitude*math.cos(s+omega*t+relative_phase)
        rows.append((s, forward, reverse, forward+reverse))
    combined = [r[3] for r in rows]
    return {
        "model": "kinematic counter-propagating helical phase, not a physical DEE evolution",
        "helical_phase": "s=m*phi+kz*z",
        "forward": "A*cos(s-omega*t)",
        "reverse": "C*cos(s+omega*t+delta)",
        "m":m, "kz":kz, "omega":omega, "time":t,
        "relative_phase":relative_phase,
        "combined_peak_abs":max(abs(v) for v in combined),
        "combined_rms":math.sqrt(sum(v*v for v in combined)/n),
        "samples": [{"s":s,"forward":f,"reverse":r,"combined":c}
                    for s,f,r,c in rows],
        "limitations": "No density, pressure, stress-energy, dispersion relation, or constraint solve; compression requires a separately defined physical field."
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--m",type=int,default=2)
    p.add_argument("--kz",type=float,default=1.0)
    p.add_argument("--omega",type=float,default=1.0)
    p.add_argument("--forward-amplitude",type=float,default=1.0)
    p.add_argument("--reverse-amplitude",type=float,default=1.0)
    p.add_argument("--relative-phase",type=float,default=0.0)
    p.add_argument("--time",type=float,default=0.0)
    p.add_argument("--samples",type=int,default=256)
    args=p.parse_args()
    print(json.dumps(sample(args.m,args.kz,args.omega,args.forward_amplitude,
                            args.reverse_amplitude,args.relative_phase,
                            args.samples,args.time),indent=2))


if __name__=="__main__":
    main()
