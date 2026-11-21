"""Decompose the frozen inner RMS without changing its norm or sample set."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'reproduce'))
from diagnose_toroidal_boundary_convergence import Audit


def analyze(directory):
    records = []
    for n in (8,16,32,64):
        with np.load(directory/f'fields_{n}.npz') as fields:
            u,w = fields['psi'],fields['W_cartesian']
        a = Audit(u.shape)
        weights = a.r*a.q*u**6
        residual = a.momentum(u,w,2)
        numerator = np.sum(weights*np.sum(residual**2,axis=0),axis=(1,2))
        mass = weights.sum(axis=(1,2))
        inner = (a.r[:,0,0]<=.10+1e-12)
        inner[0] = False
        records.append(dict(n=n,inner_row_rms=np.sqrt(numerator[inner]/mass[inner]).tolist(),
                       first_row_share=float(numerator[1]/numerator[inner].sum())))
    return records


if __name__ == '__main__':
    print(json.dumps(analyze(Path(sys.argv[1])),indent=2))
