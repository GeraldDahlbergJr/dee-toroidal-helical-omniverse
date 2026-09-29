"""Boundary-distance/refinement diagnostic for the pinned reduced ADM evolution.

No continuum equations, initial data, RK4 update, or frozen two-cell boundary
treatment are changed.  Measure Hamiltonian and momentum constraints using both
fixed cell trims and matched physical boundary distances across refinement.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import _pinned_coupled_adm_smoke as adm

CELL_TRIMS=(2,4,8,12,16,24,32,48,64)
PHYSICAL_FRACTIONS=(0.01,0.02,0.04,0.08,0.12,0.16)

def vecmag(M):
    return np.sqrt(np.sum(M*M,axis=1))

def diag(H,M,r,trim):
    sl=slice(trim,-trim)
    ah=np.abs(H[sl]); am=vecmag(M[sl])
    ih=int(np.argmax(ah))+trim; im=int(np.argmax(am))+trim
    return {
      "trim_cells":int(trim),
      "trim_distance":float(trim*(r[1]-r[0])),
      "H_max":float(ah.max()),"H_r_at_max":float(r[ih]),
      "M_max":float(am.max()),"M_r_at_max":float(r[im]),
      "M_components_at_max":[float(x) for x in M[im]],
      "H_rms":float(np.sqrt(np.mean(H[sl]**2))),
      "M_rms":float(np.sqrt(np.mean(am**2))),
    }

def physical_diags(H,M,r):
    h=float(r[1]-r[0]); width=float(r[-1]-r[0]); out=[]
    for frac in PHYSICAL_FRACTIONS:
        distance=frac*width
        trim=max(2,int(np.ceil(distance/h)))
        if 2*trim < len(r):
            d=diag(H,M,r,trim)
            d["target_fraction_of_domain"]=frac
            d["target_distance"]=distance
            out.append(d)
    return out

def run(n,t_end,cfl):
    r,state=adm.initial_state(n); initial=tuple(x.copy() for x in state)
    h=r[1]-r[0]; dt0=cfl*h; steps=int(np.ceil(t_end/dt0)); dt=t_end/steps
    H0,M0,*_=adm.constraints(*state,h)
    initial_cells=[diag(H0,M0,r,t) for t in CELL_TRIMS if 2*t<n]
    initial_phys=physical_diags(H0,M0,r)
    for _ in range(steps):
        state=adm.rk4_step(state,dt,h,initial)
    Hf,Mf,*_=adm.constraints(*state,h)
    final_cells=[diag(Hf,Mf,r,t) for t in CELL_TRIMS if 2*t<n]
    final_phys=physical_diags(Hf,Mf,r)
    return {"points":n,"h":float(h),"domain":[float(r[0]),float(r[-1])],
            "dt":float(dt),"steps":steps,
            "initial":{"cell_trims":initial_cells,"physical_trims":initial_phys},
            "final":{"cell_trims":final_cells,"physical_trims":final_phys},
            "finite":bool(all(np.all(np.isfinite(x)) for x in state))}

def convergence(rows,key="M_max"):
    out=[]
    for frac in PHYSICAL_FRACTIONS:
        vals=[]
        for row in rows:
            ds=row["final"]["physical_trims"]
            hit=next((d for d in ds if d["target_fraction_of_domain"]==frac),None)
            if hit: vals.append((row["points"],hit[key]))
        orders=[]
        for i in range(len(vals)-1):
            if vals[i+1][1]>0 and vals[i][1]>0:
                orders.append(float(np.log(vals[i][1]/vals[i+1][1])/np.log(2.0)))
        out.append({"target_fraction_of_domain":frac,"values":vals,
                    "pairwise_orders":orders})
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--resolutions",nargs="+",type=int,default=[201,401,801])
    ap.add_argument("--t-end",type=float,default=.025)
    ap.add_argument("--cfl",type=float,default=.10)
    ap.add_argument("--output",default="coupled_einstein_dee_momentum_diagnostic.json")
    a=ap.parse_args()
    rows=[run(n,a.t_end,a.cfl) for n in a.resolutions]
    out={"checkpoint_name":"Reduced ADM boundary-distance/refinement diagnostic",
         "scope":"diagnostic only; pinned equations and frozen boundary unchanged",
         "question":"Does constraint contamination shrink with grid spacing or persist at matched physical distance?",
         "t_end":a.t_end,"cfl":a.cfl,
         "cell_trims":list(CELL_TRIMS),
         "physical_fractions":list(PHYSICAL_FRACTIONS),
         "runs":rows,
         "matched_physical_convergence":{
             "M_max":convergence(rows,"M_max"),
             "H_max":convergence(rows,"H_max")}}
    Path(a.output).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__": main()
