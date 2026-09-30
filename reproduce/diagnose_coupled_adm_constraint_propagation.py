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

def snapshot(state,r,h,t):
    H,M,*_=adm.constraints(*state,h)
    return {"time":float(t),"physical_trims":physical_diags(H,M,r)}

def run(n,t_end,cfl,sample_times):
    r,state=adm.initial_state(n); initial=tuple(x.copy() for x in state)
    h=float(r[1]-r[0]); dt0=cfl*h; steps=int(np.ceil(t_end/dt0)); dt=t_end/steps
    targets=sorted(set([0.0]+[float(x) for x in sample_times if 0 < x <= t_end]+[float(t_end)]))
    snaps=[]; ti=0
    if targets and targets[0]==0.0:
        snaps.append(snapshot(state,r,h,0.0)); ti=1
    for step in range(1,steps+1):
        state=adm.rk4_step(state,dt,h,initial)
        t=step*dt
        while ti < len(targets) and t+0.5*dt >= targets[ti]:
            snaps.append(snapshot(state,r,h,t)); ti+=1
    return {"points":n,"h":h,"dt":float(dt),"steps":steps,
            "finite":bool(all(np.all(np.isfinite(x)) for x in state)),
            "snapshots":snaps}

def convergence(rows,key):
    times=sorted({round(s["time"],12) for row in rows for s in row["snapshots"]})
    out=[]
    for t in times:
        by_frac=[]
        for frac in PHYSICAL_FRACTIONS:
            vals=[]
            for row in rows:
                s=min(row["snapshots"],key=lambda x:abs(x["time"]-t))
                d=next(x for x in s["physical_trims"] if x["target_fraction_of_domain"]==frac)
                vals.append((row["points"],d[key]))
            orders=[float(np.log(vals[i][1]/vals[i+1][1])/np.log(2.0))
                    for i in range(len(vals)-1) if vals[i][1]>0 and vals[i+1][1]>0]
            by_frac.append({"target_fraction_of_domain":frac,"values":vals,"pairwise_orders":orders})
        out.append({"time":t,"physical_convergence":by_frac})
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--resolutions",nargs="+",type=int,default=[201,401,801])
    ap.add_argument("--t-end",type=float,default=.20)
    ap.add_argument("--cfl",type=float,default=.10)
    ap.add_argument("--sample-times",nargs="+",type=float,default=[.025,.05,.10,.15,.20])
    ap.add_argument("--output",default="coupled_einstein_dee_constraint_propagation.json")
    a=ap.parse_args()
    rows=[run(n,a.t_end,a.cfl,a.sample_times) for n in a.resolutions]
    out={"checkpoint_name":"Reduced ADM time-resolved constraint-propagation diagnostic",
         "scope":"diagnostic only; pinned equations, initial data, RK4/CFL and frozen boundary unchanged",
         "question":"How does boundary-associated constraint contamination propagate inward with physical time and matched physical distance?",
         "t_end":a.t_end,"cfl":a.cfl,"sample_times":a.sample_times,
         "physical_fractions":list(PHYSICAL_FRACTIONS),"runs":rows,
         "matched_physical_time_convergence":{"M_max":convergence(rows,"M_max"),
                                              "H_max":convergence(rows,"H_max")}}
    Path(a.output).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__": main()
