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
    snaps=[]; ti=0
    if targets and targets[0]==0.0:
        snaps.append(snapshot(state,r,h,0.0,0.0)); ti=1
    for step in range(1,steps+1):
        state=adm.rk4_step(state,dt,h,initial)
        t=step*dt
        while ti < len(targets) and t+0.5*dt >= targets[ti]:
            snaps.append(snapshot(state,r,h,t,targets[ti])); ti+=1
    return {"points":n,"h":h,"dt":float(dt),"steps":steps,
            "finite":bool(all(np.all(np.isfinite(x)) for x in state)),
            "snapshots":snaps}

def convergence(rows,key):
    targets=sorted({s["target_time"] for row in rows for s in row["snapshots"]})
    out=[]
    for target in targets:
        by_frac=[]
        selected=[]
        for row in rows:
            s=next(x for x in row["snapshots"] if x["target_time"]==target)
            selected.append(s)
        for frac in PHYSICAL_FRACTIONS:
            vals=[]
            actual_times=[]
            offsets=[]
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

