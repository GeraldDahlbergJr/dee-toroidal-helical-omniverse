"""Physical nodal-SBP operator diagnostics, independent of a coupled solve."""
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'research'))
from run43_mixed_flux_sbp_physical import MixedFluxSBPGrid, RMIN, RMAX


def probe(n):
    g = MixedFluxSBPGrid(n)
    k = np.pi/(RMAX-RMIN)
    f = k*k*(g.r-RMIN)*(RMAX-g.r)
    fp = k*k*(RMAX+RMIN-2*g.r)
    fpp = -2*k*k
    c,s = np.cos(g.t),np.sin(g.t)
    force = (g.basis[0]*(s*fpp) + g.basis[1]*(c*fp/g.r))/3
    force[2] += fpp+(1/g.r+c/g.q)*fp
    w = np.zeros((3,g.size))
    w[2] = f
    wb = w[:,g.idx].ravel()
    err = (g.vector@wb).reshape(3,-1)-force[:,g.idx]
    rng = np.random.default_rng(43)
    u,v = rng.normal(size=(2,wb.size))
    mass = np.tile(g.sbp_weight[g.idx],3)
    uv = u@(mass*(g.vector@v))
    vu = v@(mass*(g.vector@u))
    # Weighted weak flux against the strain pairing, including angular damping separately.
    no_damp = MixedFluxSBPGrid(n,stabilization=0)
    def strain(x):
        full = np.zeros_like(w)
        full[:,g.idx] = x.reshape(3,-1)
        return no_damp.longitudinal(full)
    weak = u@(mass*(no_damp.vector@v))
    expected = -.5*np.sum(g.sbp_weight*strain(u)*strain(v))
    row = np.sqrt((err**2).sum(axis=0).reshape(n-1,n,n).mean(axis=(1,2)))
    return dict(n=n,sbp=g.radial_sbp_defect,adjoint_relative=abs(weak-expected)/abs(expected),
                stabilized_skew_absolute=abs(uv-vu),manufactured_first_rms=row[0],
                manufactured_last_rms=row[-1],manufactured_bulk_rms=row[n//2-1],
                tested_energy=float(-u@(mass*(g.vector@u))))


if __name__=='__main__':
    os.environ['DEE_PRECONDITIONER']='fourier'
    os.environ['DEE_SCALAR_ORDER']='flux'
    records=[probe(n) for n in (8,16,32,64)]
    Path(sys.argv[1]).write_text(json.dumps(records,indent=2)+'\n')
