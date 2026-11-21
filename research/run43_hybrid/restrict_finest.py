"""Sample one finest physical solution on coarser audit nodes to isolate field changes."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'reproduce'))
from diagnose_toroidal_boundary_convergence import Audit


def compare(directory):
    with np.load(directory/'fields_64.npz') as fields:
        u64,w64 = fields['psi'],fields['W_cartesian']
    rows = []
    for n in (8,16,32):
        with np.load(directory/f'fields_{n}.npz') as fields:
            u,w = fields['psi'],fields['W_cartesian']
        a = Audit(u.shape)
        step = 64//n
        coarse = a.momentum(u,w,2)
        restricted = a.momentum(u64[::step,::step,::step],w64[:,::step,::step,::step],2)
        mask = np.broadcast_to((a.r>a.r[0])&(a.r<=.10+1e-12),u.shape)
        def norm(v):
            return a.norms(v,u,mask)['volume_weighted_rms']
        rows.append(dict(n=n,coarse_rms=norm(coarse),restricted_rms=norm(restricted),
                         difference_rms=norm(coarse-restricted)))
    return rows


if __name__=='__main__':
    print(json.dumps(compare(Path(sys.argv[1])),indent=2))
