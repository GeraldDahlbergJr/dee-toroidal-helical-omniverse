"""Manufactured inner-boundary check of the independent momentum audit.

Tests the derivative operator, not the Run43 candidate or its physics.
For W=(x^2,y^2,z^2), div(L W)=(8/3,8/3,8/3) in Cartesian space.
No solver/checkpoint modification. Run from reproduce/ or repository root.
"""
import json
import numpy as np
from diagnose_toroidal_boundary_convergence import Audit
from validate_toroidal_helical_dee_source import R, RMIN

def check(shape):
    audit = Audit(shape)
    x = audit.q * np.cos(audit.p) + np.zeros(shape)
    y = audit.q * np.sin(audit.p) + np.zeros(shape)
    z = audit.r * np.sin(audit.t) + np.zeros(shape)
    w = np.array([x*x, y*y, z*z])
    expected = np.full((3,) + shape, 8.0/3.0)
    mask = np.broadcast_to(((audit.r[:,0,0] >= RMIN) &
                            (audit.r[:,0,0] <= .10))[:,None,None], shape).copy()
    mask[0] = False
    results = {}
    for order in (2,4):
        a = audit.tensor(w, order)
        divergence = np.array([
            sum(audit.grad(a[i,j], order)[j] for j in range(3))
            for i in range(3)
        ])
        diff = divergence - expected
        results[str(order)] = {
            'inner_band_rms': float(np.sqrt(np.mean(np.sum(diff[:,mask]**2,axis=0)))),
            'inner_band_max': float(np.max(np.linalg.norm(diff[:,mask],axis=0))),
            'interior_samples': int(mask.sum()),
        }
    return {'shape': list(shape), 'errors': results}

def main():
    rows = [check((n, n-1, n-1)) for n in (9,17,33,65)]
    observed = {}
    for order in ('2','4'):
        e = [r['errors'][order]['inner_band_rms'] for r in rows]
        observed[order] = [float(np.log2(e[i]/e[i+1])) for i in range(len(e)-1)]
    report = {
        'purpose':'Manufactured polynomial check of independent momentum differential operator near inner boundary',
        'exact_divergence_LW':[8/3]*3,
        'rows':rows,'observed_orders':observed,
        'caution':'A successful manufactured operator test does not validate Run43 physical boundary conditions.'
    }
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
