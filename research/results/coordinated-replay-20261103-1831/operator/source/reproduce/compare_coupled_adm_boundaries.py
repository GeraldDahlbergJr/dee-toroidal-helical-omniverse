"""Controlled boundary-treatment comparison for the pinned reduced ADM evolution.

No continuum equations are changed. Compare the legacy frozen two-cell buffer
against simple second-order zero-normal-gradient extrapolation of evolved
variables after each RK4 step. This is a diagnostic, not a production boundary
condition and not a claim of constraint preservation.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import _pinned_coupled_adm_smoke as adm

def extrapolate(u):
    v=u.copy()
    # linear extrapolation from interior: second difference zero at boundary cells
    v[1]=2*v[2]-v[3]; v[0]=2*v[1]-v[2]
    v[-2]=2*v[-3]-v[-4]; v[-1]=2*v[-2]-v[-3]
    return v

def step_extrap(state,dt,h,initial):
    k1=adm.rhs(state,h)
    k2=adm.rhs(adm.add_state(state,k1,.5*dt),h)
    k3=adm.rhs(adm.add_state(state,k2,.5*dt),h)
    k4=adm.rhs(adm.add_state(state,k3,dt),h)
    out=tuple(x+dt*(a+2*b+2*c+d)/6 for x,a,b,c,d in zip(state,k1,k2,k3,k4))
    return tuple(extrapolate(x) for x in out)

def norms_components(H,M,trim):
    sl=slice(trim,-trim); mag=np.sqrt(np.sum(M[sl]**2,axis=1))
    return {"H_max":float(np.max(np.abs(H[sl]))),
      "Mr_max":float(np.max(np.abs(M[sl,0]))),
      "Mphi_max":float(np.max(np.abs(M[sl,1]))),
      "Mz_max":float(np.max(np.abs(M[sl,2]))),
      "Mmag_max":float(np.max(mag))}

def run(n,t_end,cfl,mode):
    r,state=adm.initial_state(n); initial=tuple(x.copy() for x in state)
    h=r[1]-r[0]; steps=int(np.ceil(t_end/(cfl*h))); dt=t_end/steps
    H0,M0,*_=adm.constraints(*state,h)
    for _ in range(steps):
        state=(adm.rk4_step(state,dt,h,initial) if mode=="frozen2"
               else step_extrap(state,dt,h,initial))
    H,M,*_=adm.constraints(*state,h)
    return {"points":n,"dt":dt,"steps":steps,
      "initial":{str(t):norms_components(H0,M0,t) for t in (8,16,32)},
      "final":{str(t):norms_components(H,M,t) for t in (8,16,32)},
      "finite":bool(all(np.all(np.isfinite(x)) for x in state))}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--resolutions",nargs="+",type=int,default=[201,401,801])
    ap.add_argument("--t-end",type=float,default=.025); ap.add_argument("--cfl",type=float,default=.1)
    ap.add_argument("--output",default="coupled_einstein_dee_boundary_comparison.json"); a=ap.parse_args()
    out={"checkpoint_name":"Reduced ADM boundary-treatment comparison",
      "scope":"diagnostic only; same continuum equations and initial data",
      "warning":"extrapolated treatment is a numerical comparison, not asserted to be a mathematically constraint-preserving boundary condition",
      "t_end":a.t_end,"cfl":a.cfl,"treatments":{}}
    for mode in ("frozen2","extrapolated2"):
        out["treatments"][mode]=[run(n,a.t_end,a.cfl,mode) for n in a.resolutions]
    Path(a.output).write_text(json.dumps(out,indent=2)+"\n"); print(json.dumps(out,indent=2))
if __name__=="__main__": main()
