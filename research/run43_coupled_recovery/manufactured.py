"""Analytic tensor forcing and a discrete wall-lifting regression."""
import argparse
import json
from pathlib import Path

import numpy as np

from solver import GRIDS, Momentum, RMIN, RMAX
from diagnose import Audit, R, basis, corrections, frozen_parts
from operator_eval import divergence


def measure(n):
    op = Momentum(n)
    x = (op.r-RMIN)/(RMAX-RMIN)
    theta = op.theta
    cases = {}
    for name, phase, frequency in [('homogeneous', 0., 1), ('lifted', .4, 2)]:
        k = np.pi/(RMAX-RMIN)
        f = np.sin(np.pi*x+phase)
        if phase == 0.:
            f[[0, -1]] = 0.
        fp = k*np.cos(np.pi*x+phase)
        g = np.concatenate([np.cos(frequency*theta), np.sin(theta), 1+.2*np.cos(theta)])
        gt = np.concatenate([-frequency*np.sin(frequency*theta), np.cos(theta), -.2*np.sin(theta)])
        gtt = np.concatenate([-frequency**2*np.cos(frequency*theta), -np.sin(theta), -.2*np.cos(theta)])
        field = f*g[:, None]
        exact = divergence(field, fp*g[:, None], -k*k*field, f*gt[:, None],
                           fp*gt[:, None], f*gtt[:, None], op.r, theta, R)
        actual = op.apply(field)
        solved, info = op.solve(exact, field)
        discrete, discrete_info = op.solve(actual, field)
        np.testing.assert_allclose(discrete, field, atol=2e-12, rtol=0)
        np.testing.assert_array_equal(solved[:, [0, -1]], field[:, [0, -1]])

        # Verify all ordered terms against the prior substitution on free rows.
        a, da = basis(op.r, theta)
        audit = Audit((n+1, n, n))
        parts = frozen_parts(audit, field, a)
        expected = sum(parts.values())+corrections(audit, field, a, da, parts)['compact_rr_and_metric']
        np.testing.assert_allclose(actual[:, 1:-1], expected[:, 1:-1], atol=1e-9, rtol=1e-12)

        def rms(v):
            return float(np.sqrt(np.mean(np.sum(v*v, axis=0))))

        regions = dict(inner_wall=0, inner_first=1, interior=slice(2, -2),
                       outer_first=-2, outer_wall=-1)
        cases[name] = dict(
            truncation={label: rms((actual-exact)[:, region]) for label, region in regions.items()},
            field_error={label: rms((solved-field)[:, region]) for label, region in regions.items()},
            discrete_recovery_linf=float(np.max(abs(discrete-field))),
            continuous_solve=info, discrete_solve=discrete_info,
        )
    return dict(n=n, cases=cases, stability=op.stability())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    rows = [measure(n) for n in (*GRIDS, 128)]
    for case in rows[-1]['cases']:
        for region, error in rows[-1]['cases'][case]['truncation'].items():
            prior = rows[-2]['cases'][case]['truncation'][region]
            assert np.log2(prior/error) > 1.8, (case, region, prior, error)
    args.output.write_text(json.dumps(rows, indent=2)+'\n')


if __name__ == '__main__':
    main()
