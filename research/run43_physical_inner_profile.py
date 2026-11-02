"""Independent physical-field radial residual profiling; diagnostic only.

Uses original repository n=8,16,32 coupled solutions without modifying their
sources, boundary conditions, frozen audit or acceptance thresholds.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "reproduce"))
from diagnose_toroidal_boundary_convergence import Audit
from validate_toroidal_helical_dee_source import RMIN, RMAX

def main():
    rows = []
    for n in (8, 16, 32):
        path = ROOT / "reproduce" / f"toroidal_coupled_fields_{n+1}_{n}.npz"
        with np.load(path) as data:
            u = data["psi"]
            w = data["W_cartesian"]
        assert u.shape == (n+1,n,n) and w.shape == (3,n+1,n,n)
        assert np.isfinite(u).all() and np.isfinite(w).all() and np.min(u)>0
        audit = Audit(u.shape)
        h = audit.flux_hamiltonian(u,w)
        m = audit.momentum(u,w,2)
        weight = np.broadcast_to(audit.r*audit.q,u.shape)*u**6
        profile=[]
        for i in range(1,n):
            sl=np.s_[i,:,:]
            hn=float(np.sqrt(np.sum(h[sl]**2*weight[sl])/np.sum(weight[sl])))
            mn=float(np.sqrt(np.sum(np.sum(m[:,i]**2,axis=0)*weight[sl])/np.sum(weight[sl])))
            profile.append({"r":float(audit.r[i,0,0]),"H_rms":hn,"M2_rms":mn})
        rows.append({"n":n,"profile":profile,
                     "inner":audit.summarize(u,w)["inner_boundary_band"]})
    result={"scope":"Original coupled physical coarse fields only; no n=64 result or 12/12 validation",
            "radial_profiles":rows}
    print(json.dumps(result,indent=2))
    for a,b in zip(rows,rows[1:]):
        for metric in ("independent_flux_H","independent_tensor_M_order2"):
            x=a["inner"][metric]["volume_weighted_rms"]
            y=b["inner"][metric]["volume_weighted_rms"]
            print(f"{metric} {a['n']}->{b['n']}: p={np.log2(x/y):.6f} (original physical coarse fields)")
    print("Diagnostic only: do not infer the 32->64 gate or modify frozen acceptance.")

if __name__=="__main__":
    main()
