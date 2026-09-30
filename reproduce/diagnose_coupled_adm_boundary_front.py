"""Diagnostic-only boundary-contamination-front measurement for the pinned ADM gate."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import _pinned_coupled_adm_smoke as adm
FRACTIONS=(.01,.02,.03,.04,.05,.06,.07,.08,.10,.12,.14,.16,.18,.20)
ORDER_FLOOR=1.7
def mag(M): return np.sqrt(np.sum(M*M,axis=1))
def diags(H,M,r):
 h=float(r[1]-r[0]); w=float(r[-1]-r[0]); out=[]
 for f in FRACTIONS:
  trim=max(2,int(np.ceil(f*w/h))); sl=slice(trim,-trim)
  out.append({"fraction":f,"trim_cells":trim,"trim_distance":trim*h,
              "H_max":float(np.max(np.abs(H[sl]))),"M_max":float(np.max(mag(M[sl])))})
 return out
def snap(state,r,h,t,target):
 H,M,*_=adm.constraints(*state,h)
 return {"target_time":float(target),"actual_time":float(t),"time_offset":float(t-target),
         "physical_trims":diags(H,M,r)}
def run(n,t_end,cfl,times):
 r,state=adm.initial_state(n); initial=tuple(x.copy() for x in state); h=float(r[1]-r[0])
 steps=int(np.ceil(t_end/(cfl*h))); dt=t_end/steps
 targets=sorted(set([0.]+[float(x) for x in times if 0<x<=t_end]+[float(t_end)]))
 snaps=[snap(state,r,h,0.,0.)]; ti=1; t=0.
 while t<t_end:
  target=targets[ti] if ti<len(targets) else t_end; d=min(dt,target-t,t_end-t)
  if d<=0:
   if ti<len(targets): t=targets[ti]; snaps.append(snap(state,r,h,t,targets[ti])); ti+=1; continue
   break
  state=adm.rk4_step(state,d,h,initial); t+=d
  tol=32*np.finfo(float).eps*max(1.,abs(t))
  if ti<len(targets) and abs(t-targets[ti])<=tol:
   t=targets[ti]; snaps.append(snap(state,r,h,t,targets[ti])); ti+=1
 return {"points":n,"h":h,"dt":dt,"steps":steps,
         "finite":bool(all(np.all(np.isfinite(x)) for x in state)),"snapshots":snaps}
def order(a,b): return None if a<=0 or b<=0 else float(np.log(a/b)/np.log(2.))
def analyze(rows):
 by={r["points"]:r for r in rows}; out=[]
 for a in by[401]["snapshots"]:
  t=a["target_time"]; b=next(s for s in by[801]["snapshots"] if s["target_time"]==t); fs=[]
  for x in a["physical_trims"]:
   y=next(d for d in b["physical_trims"] if d["fraction"]==x["fraction"])
   pH=order(x["H_max"],y["H_max"]); pM=order(x["M_max"],y["M_max"])
   fs.append({"fraction":x["fraction"],"p_H":pH,"p_M":pM,
              "clean_second_order":bool(pH is not None and pM is not None and pH>=ORDER_FLOOR and pM>=ORDER_FLOOR)})
  clean=[x["fraction"] for x in fs if x["clean_second_order"]]
  out.append({"target_time":t,"innermost_clean_fraction":min(clean) if clean else None,
              "fine_pair_order_floor":ORDER_FLOOR,"fractions":fs})
 return out
def main():
 p=argparse.ArgumentParser(); p.add_argument("--resolutions",nargs="+",type=int,default=[201,401,801])
 p.add_argument("--t-end",type=float,default=.4); p.add_argument("--cfl",type=float,default=.1)
 p.add_argument("--sample-times",nargs="+",type=float,default=[.05,.10,.15,.20,.25,.30,.35,.40])
 p.add_argument("--output",default="coupled_einstein_dee_boundary_front.json"); a=p.parse_args()
 rows=[run(n,a.t_end,a.cfl,a.sample_times) for n in a.resolutions]
 finite=all(r["finite"] for r in rows); exact=all(s["time_offset"]==0. for r in rows for s in r["snapshots"])
 payload={"checkpoint_name":"Coupled Einstein-DEE boundary-contamination-front gate",
 "status":"MEASURED" if finite and exact else "INCONCLUSIVE",
 "scope":"diagnostic only; pinned reduced cylindrical-annulus ADM equations, initial data, RK4/CFL ceiling, refinement sequence, and frozen two-cell boundary treatment unchanged",
 "t_end":a.t_end,"cfl":a.cfl,"resolutions":a.resolutions,"sample_times":a.sample_times,
 "fractions":list(FRACTIONS),"order_floor":ORDER_FLOOR,"all_runs_finite":finite,
 "exact_time_alignment":exact,"runs":rows,"front_analysis":analyze(rows),
 "interpretation":"Innermost matched physical trim where both Hamiltonian and momentum 401->801 orders are >=1.7; diagnostic only, not a causal-speed or physical-boundary claim."}
 Path(a.output).write_text(json.dumps(payload,indent=2)+"\n"); print(json.dumps(payload,indent=2))
 if not(finite and exact): raise SystemExit(3)
if __name__=="__main__": main()
