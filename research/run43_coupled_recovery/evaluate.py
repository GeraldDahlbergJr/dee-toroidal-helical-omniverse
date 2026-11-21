"""Evaluate frozen, substituted, and continuum momentum on saved fields."""
import argparse
import csv
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from solver import GRIDS, RMIN, RMAX
from diagnose import (ARCHIVE, MANIFEST, Audit, R, G, basis, corrections,
                      digest, frozen_parts, reference_fields, restrictions,
                      rms)
from operator_eval import traces

HYPOTHESIS = 'compact_rr_and_metric'


def quartic_samples(w, bias):
    """Resample each coarse interval from a local quartic for the trace evaluator."""
    n = w.shape[1]-1
    result = np.empty((3, 8*n+1, n))
    points = np.arange(5)
    for i in range(n):
        start = min(max(i-bias, 0), n-4)
        x = i-start+np.arange(1, 8)/8
        coefficients = np.column_stack([
            np.prod((x[:, None]-np.delete(points, j))/(j-np.delete(points, j)), axis=1)
            for j in points
        ])
        result[:, 8*i+1:8*(i+1)] = np.einsum('ij,cjt->cit', coefficients, w[:, start:start+5])
    result[:, ::8] = w
    return result


def load(directory, n):
    raw = (directory/f'fields_{n}.npz').read_bytes()
    expected = json.loads((directory/f'result_{n}.json').read_text())['sha256']
    assert digest(raw) == expected
    with np.load(directory/f'fields_{n}.npz') as data:
        return data['psi'][..., 0], data['W_cartesian'][..., 0]


def weighted_dot(weight, a, b):
    return np.sum(weight*np.sum(a*b, axis=0))/weight.sum()


def measure(directory, prior):
    rows = []
    frozen_report = json.loads((directory/'independent-audit.json').read_text())
    with ZipFile(ARCHIVE) as archive:
        for n, frozen_row in zip(GRIDS, frozen_report['rows']):
            hybrid = load(prior/'hybrid', n)
            assert digest((prior/'hybrid'/f'fields_{n}.npz').read_bytes()) == MANIFEST['hybrid_sha256'][str(n)]
            uref, wref = reference_fields(archive, prior/'levels', n, 4, {})
            audit = Audit((n+1, n, n))
            r, theta = audit.r[..., 0], audit.t[..., 0]
            a, da = basis(r, theta)
            for method, (u, w) in [('coupled', load(directory, n)), ('hybrid', hybrid),
                                    ('reference_4N', (uref, wref[:, ::8, ::4]))]:
                parts = frozen_parts(audit, w, a)
                fixes = corrections(audit, w, a, da, parts)
                source = 8*np.pi*G*u**-4*audit.S[..., 0]
                frozen = u**-10*sum(parts.values())-source
                combined = frozen+u**-10*fixes[HYPOTHESIS]
                rr_gap = -u**-10*fixes['compact_rr']
                samples = wref if method == 'reference_4N' else quartic_samples(w, 1)
                alternative = samples if method == 'reference_4N' else quartic_samples(w, 2)
                regions = {**restrictions(audit), "interior": np.arange(1, n)}
                for region, ix in regions.items():
                    weight = (r*audit.q[..., 0]*u**6)[ix]
                    def continuum(field):
                        div = sum(traces(field, n, ix, r[ix], theta, R, RMAX-RMIN))/2
                        return u[ix]**-10*div-source[:, ix]

                    cont = continuum(samples)
                    cont_alternative = continuum(alternative)
                    gap = frozen[:, ix]-cont
                    diff = w[:, ix]-wref[:, ::8, ::4][:, ix]
                    du = u[ix]-uref[ix]
                    row = dict(n=n, field=method, region=region,
                               frozen_rms=rms(weight, frozen[:, ix]),
                               combined_rms=rms(weight, combined[:, ix]),
                               continuum_rms=rms(weight, cont),
                               alternate_continuum_rms=rms(weight, cont_alternative),
                               reconstruction_delta_rms=rms(weight, cont_alternative-cont),
                               frozen_continuum_gap_rms=rms(weight, gap),
                               rr_gap_rms=rms(weight, rr_gap[:, ix]),
                               gap_after_rr_rms=rms(weight, gap-rr_gap[:, ix]),
                               rr_gap_projection=weighted_dot(weight, rr_gap[:, ix], gap)/weighted_dot(weight, gap, gap),
                               W_reference_rms=rms(weight, diff),
                               psi_reference_rms=rms(weight, du[None]),
                               W_hybrid_rms=rms(weight, w[:, ix]-hybrid[1][:, ix]),
                               psi_hybrid_rms=rms(weight, (u[ix]-hybrid[0][ix])[None]))
                    if method == 'coupled' and region.endswith('_band'):
                        key = region.replace('_band', '_boundary_band')
                        expected = frozen_row['audit'][key]['independent_tensor_M_order2']['volume_weighted_rms']
                        np.testing.assert_allclose(row['frozen_rms'], expected, atol=1e-15, rtol=0)
                        assert row['combined_rms'] < 1e-13
                    rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('prior', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    rows = measure(args.directory, args.prior)
    previous = {}
    for row in rows:
        key = (row['field'], row['region'])
        for label in ('frozen', 'combined', 'continuum', 'W_reference', 'psi_reference'):
            value = row[f'{label}_rms']
            old = previous.get(key, {}).get(f'{label}_rms', 0)
            row[f'{label}_order'] = float(np.log2(old/value)) if old > 0 and value > 0 else ''
        previous[key] = row
    with args.output.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
