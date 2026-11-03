"""Time-resolved constraint-propagation diagnostic for the pinned reduced ADM evolution.

No continuum equations, initial data, RK4 update, CFL, refinement sequence, or
frozen two-cell boundary treatment are changed. This diagnostic samples the same
evolution at fixed physical times and matched physical boundary distances.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import _pinned_coupled_adm_smoke as adm
from diagnose_coupled_adm_momentum import PHYSICAL_FRACTIONS, physical_diags

def snapshot(state,r,h,t,target_time):
    H,M,*_=adm.constraints(*state,h)
    return {"target_time":float(target_time),"actual_time":float(t),
            "time_offset":float(t-target_time),
            "physical_trims":physical_diags(H,M,r)}

def run(n,t_end,cfl,sample_times):
    r,state=adm.initial_state(n); initial=tuple(x.copy() for x in state)
    h=float(r[1]-r[0]); dt0=cfl*h; steps=int(np.ceil(t_end/dt0)); dt=t_end/steps
    targets=sorted(set([0.0]+[float(x) for x in sample_times if 0 < x <= t_end]+[float(t_end)]))
    snaps=[]; ti=0; t=0.0
    if targets and targets[0]==0.0:
        snaps.append(snapshot(state,r,h,0.0,0.0)); ti=1
    # Keep the nominal CFL step as a ceiling, but shorten the step whenever
    # necessary so every requested diagnostic time is hit exactly.
    while t < t_end:
        next_target=targets[ti] if ti < len(targets) else t_end
        step_dt=min(dt, next_target-t, t_end-t)
        if step_dt <= 0.0:
            if ti < len(targets) and abs(t-targets[ti]) <= 32*np.finfo(float).eps*max(1.0,abs(t)):
                snaps.append(snapshot(state,r,h,targets[ti],targets[ti])); ti+=1
                continue
            break
        state=adm.rk4_step(state,step_dt,h,initial)
        t += step_dt
        tol=32*np.finfo(float).eps*max(1.0,abs(t))
        if ti < len(targets) and abs(t-targets[ti]) <= tol:
            t=targets[ti]
            snaps.append(snapshot(state,r,h,t,targets[ti])); ti+=1
    return {"points":n,"h":h,"dt":float(dt),"steps":steps,
            "finite":bool(all(np.all(np.isfinite(x)) for x in state)),
            "max_nominal_dt":float(dt),
            "snapshots":snaps}

def convergence(rows,key):
    targets=sorted({s["target_time"] for row in rows for s in row["snapshots"]})
    out=[]
    for target in targets:
        by_frac=[]; selected=[]
        for row in rows:
            s=next(x for x in row["snapshots"] if x["target_time"]==target)
            selected.append(s)
        for frac in PHYSICAL_FRACTIONS:
            vals=[]; actual_times=[]; offsets=[]
            for row,s in zip(rows,selected):
                d=next(x for x in s["physical_trims"] if x["target_fraction_of_domain"]==frac)
                vals.append((row["points"],d[key]))
                actual_times.append((row["points"],s["actual_time"]))
                offsets.append((row["points"],s["time_offset"]))
            orders=[float(np.log(vals[i][1]/vals[i+1][1])/np.log(2.0))
                    for i in range(len(vals)-1) if vals[i][1]>0 and vals[i+1][1]>0]
            by_frac.append({"target_fraction_of_domain":frac,"values":vals,
                            "actual_times":actual_times,"time_offsets":offsets,
                            "pairwise_orders":orders})
        out.append({"target_time":target,"physical_convergence":by_frac})
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--resolutions",nargs="+",type=int,default=[201,401,801])
    p.add_argument("--t-end",type=float,default=0.20)
    p.add_argument("--cfl",type=float,default=0.10)
    p.add_argument("--sample-times",nargs="+",type=float,default=[0.025,0.05,0.10,0.15,0.20])
    p.add_argument("--output",default="coupled_einstein_dee_constraint_propagation.json")
    a=p.parse_args()
    rows=[run(n,a.t_end,a.cfl,a.sample_times) for n in a.resolutions]
    payload={
        "diagnostic":"time-resolved constraint propagation at matched physical boundary distances",
        "scope":"pinned reduced cylindrical-annulus ADM evolution; frozen two-cell boundary treatment",
        "t_end":a.t_end,"cfl":a.cfl,"resolutions":a.resolutions,
        "sample_times":a.sample_times,
        "runs":rows,
        "momentum_convergence":convergence(rows,"M_max"),
        "hamiltonian_convergence":convergence(rows,"H_max"),
    }
    out=Path(a.output)
    out.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps(payload,indent=2))
    print(f"Wrote {out.resolve()}")

if __name__=="__main__":
    main()
