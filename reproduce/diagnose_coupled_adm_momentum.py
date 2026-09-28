"""Diagnose momentum-constraint growth in the pinned reduced ADM evolution.

Requires reproduce/_pinned_coupled_adm_smoke.py to have been materialized by the
workflow.  No equations are changed.  This only measures where the constraint
maximum lives and how sensitive it is to progressively wider boundary trims.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import _pinned_coupled_adm_smoke as adm

TRIMS=(2,4,8,12,16,24,32)

def vecmag(M):
    return np.sqrt(np.sum(M*M,axis=1))

def diag(H,M,r,trim):
    sl=slice(trim,-trim)
    ah=np.abs(H[sl]); am=vecmag(M[sl])
    ih=int(np.argmax(ah))+trim; im=int(np.argmax(am))+trim
    return {
      "trim":trim,
      "H_max":float(ah.max()),"H_r_at_max":float(r[ih]),
      "M_max":float(am.max()),"M_r_at_max":float(r[im]),
      "M_components_at_max":[float(x) for x in M[im]],
      "H_rms":float(np.sqrt(np.mean(H[sl]**2))),
      "M_rms":float(np.sqrt(np.mean(am**2))),
    }

def run(n,t_end,cfl):
    r,state=adm.initial_state(n); initial=tuple(x.copy() for x in state)
    h=r[1]-r[0]; dt0=cfl*h; steps=int(np.ceil(t_end/dt0)); dt=t_end/steps
    H0,M0,*_=adm.constraints(*state,h)
    initial_diag=[diag(H0,M0,r,t) for t in TRIMS if 2*t<n]
    for _ in range(steps):
        state=adm.rk4_step(state,dt,h,initial)
    Hf,Mf,*_=adm.constraints(*state,h)
    final_diag=[diag(Hf,Mf,r,t) for t in TRIMS if 2*t<n]
    return {"points":n,"h":float(h),"dt":float(dt),"steps":steps,
            "initial":initial_diag,"final":final_diag}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--resolutions",nargs="+",type=int,default=[201,401,801])
    ap.add_argument("--t-end",type=float,default=.025)
    ap.add_argument("--cfl",type=float,default=.10)
    ap.add_argument("--output",default="coupled_einstein_dee_momentum_diagnostic.json")
    a=ap.parse_args()
    rows=[run(n,a.t_end,a.cfl) for n in a.resolutions]
    out={"checkpoint_name":"Reduced ADM momentum-growth localization diagnostic",
         "scope":"diagnostic only; pinned equations unchanged",
         "t_end":a.t_end,"cfl":a.cfl,"trims":list(TRIMS),"runs":rows}
    Path(a.output).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__": main()
