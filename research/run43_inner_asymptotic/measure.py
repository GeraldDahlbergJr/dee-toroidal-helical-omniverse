"""Measure frozen momentum defects without changing the acceptance evaluator."""
import argparse
from functools import partial
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from zipfile import ZipFile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for relative in ('run43_coupled_recovery', 'run43_frozen_momentum'):
    sys.path.insert(0, str(ROOT / 'research' / relative))
from evaluate import load, quartic_samples
from diagnose import (ARCHIVE, AXES, G, GRIDS, MANIFEST, PAIRS, R, RMAX, RMIN,
                      Audit, basis, continuum_parts, contract, digest, dot,
                      frozen_parts, longitude, multiply, radial_matrix,
                      reference_fields, restrictions, rms, tensor, write_table, J)


def exact_parts(r, theta, first, second):
    a, da = basis(r, theta)
    parts = {}
    for i, j in PAIRS:
        key = AXES[i] + AXES[j]
        parts[key] = contract(tensor(second[i][j], a[j]), a[i])
        parts[key+'_metric'] = contract(tensor(first[j], da[i][j]), a[i])
    return parts


def decomposition(frozen, continuum, u, source):
    defect = {key: u**-10*(frozen[key]-continuum[key]) for key in frozen}
    mixed = ('rt_metric', 'tr_metric')
    pieces = dict(rr=defect['rr'], rt_metric=defect[mixed[0]],
                  tr_metric=defect[mixed[1]],
                  other_metric=sum(v for k, v in defect.items()
                                   if k.endswith('_metric') and k not in mixed),
                  other_principal=sum(v for k, v in defect.items()
                                      if len(k) == 2 and k != 'rr'),
                  continuum=u**-10*sum(continuum.values())-source)
    actual = u**-10*sum(frozen.values())-source
    np.testing.assert_allclose(sum(pieces.values()), actual, atol=3e-11, rtol=1e-9)
    return pieces, actual


def record(rows, interactions, label, n, h, inner, ix, weights, pieces, actual, wall=False):
    groups = dict(first=np.array([0]), remaining=np.arange(1, len(inner)),
                  inner=np.arange(len(inner)))
    if wall:
        groups['wall'] = np.array([len(ix)-1])
    mass = weights[:len(inner)].sum()
    for region, positions in groups.items():
        if not len(positions):
            continue
        weight = weights[positions]
        values = {k: v[:, positions] for k, v in pieces.items()}
        residual = actual[:, positions]
        row = dict(field=label, n=n, h=h, region=region, radial_rows=len(positions),
                   mass_fraction=weight.sum()/mass,
                   frozen_rms=rms(weight, residual))
        for key, value in values.items():
            row[key+'_rms'] = rms(weight, value)
        metric = values['rt_metric']+values['tr_metric']+values['other_metric']
        row['metric_rms'] = rms(weight, metric)
        row['radial_poloidal_metric_rms'] = rms(weight, values['rt_metric']+values['tr_metric'])
        row['gap_rms'] = rms(weight, residual-values['continuum'])
        row['frozen_energy_band_normalized'] = row['mass_fraction']*row['frozen_rms']**2
        rows.append(row)
        for (a, va), (b, vb) in combinations_with_replacement(values.items(), 2):
            product = dot(weight, va, vb)
            interactions.append(dict(field=label, n=n, region=region, a=a, b=b,
                                     signed_dot=product, energy_band_normalized=
                                     row['mass_fraction']*product*(1 if a == b else 2)))


def manufactured(rows, interactions):
    for n in (*GRIDS, 128, 256, 512, 1024, 2048):
        h, angle = (RMAX-RMIN)/n, 2*np.pi/n
        m = n//4
        # Two extra rows give the last band row its unchanged composed stencil.
        # Keep the entire shell on grids used for the direct 3D checks.
        nr = n+1 if n <= 16 else m+3
        r = np.linspace(RMIN, RMAX, n+1)[:nr, None]
        theta = (angle*np.arange(n))[None, :]
        # Equivariant reduction uses the original first-derivative routines.
        reduced = SimpleNamespace(drad={2: radial_matrix(nr, h)}, dt=angle, dp=angle,
                                  angular_derivative=partial(Audit.angular_derivative, None))
        inner = np.arange(1, m+1)
        ix = np.r_[inner, 0]
        length = RMAX-RMIN
        x = (r-RMIN)/length
        g = np.concatenate((np.cos(2*theta), np.sin(theta), 1+.2*np.cos(theta)), axis=0)[:, None]
        gt = np.concatenate((-2*np.sin(2*theta), np.cos(theta), -.2*np.sin(theta)), axis=0)[:, None]
        gtt = np.concatenate((-4*np.cos(2*theta), -np.sin(theta), -.2*np.cos(theta)), axis=0)[:, None]
        u = 1+.05*x*(1-x)*np.cos(theta)
        for label, power in (('cubic', 3), ('wall_flat_third', 4)):
            f = x-x**power
            fp = (1-power*x**(power-1))/length
            fpp = -power*(power-1)*x**(power-2)/length**2
            w = f*g
            w[:, 0] = 0
            if nr == n+1:
                w[:, -1] = 0
            wr, wt = fp[ix]*g, f[ix]*gt
            wp = multiply(J, w[:, ix])
            wrt, wrp, wtp = fp[ix]*gt, multiply(J, wr), multiply(J, wt)
            exact = exact_parts(r[ix], theta, (wr, wt, wp),
                                ((fpp[ix]*g, wrt, wrp), (wrt, f[ix]*gtt, wtp),
                                 (wrp, wtp, multiply(J, wp))))
            frozen = {k: v[:, ix] for k, v in frozen_parts(reduced, w, basis(r, theta)[0]).items()}
            source = u[ix]**-10*sum(exact.values())
            pieces, actual = decomposition(frozen, exact, u[ix], source)
            # The row-1 Taylor identity is exact on these polynomials.
            f3 = -power*(power-1)*(power-2)*x[1]**(power-3)/length**3
            f4 = -24/length**4 if power == 4 else 0
            er = basis(r[[1]], theta)[0][0]
            predicted = u[[1]]**-10*contract(tensor((h*f3/4+5*h**2*f4/24)*g, er), er)
            np.testing.assert_allclose(pieces['rr'][:, :1], predicted, atol=2e-8, rtol=1e-7)
            if n <= 16:
                audit = Audit((n+1, n, n))
                u3 = np.broadcast_to(u[..., None], audit.shape)
                w3 = longitude(w, audit.p[0, 0])
                raw = audit.momentum(u3, w3, 2)+8*np.pi*G*u3**-4*audit.S
                np.testing.assert_allclose(raw[..., 0][:, ix]-source, actual, atol=3e-11)
                raw[:, ix] -= longitude(source, audit.p[0, 0])
                mask = np.zeros(audit.shape, dtype=bool)
                mask[inner] = True
                np.testing.assert_allclose(audit.norms(raw, u3, mask)['volume_weighted_rms'],
                                           rms((r*(R+r*np.cos(theta))*u**6)[inner],
                                               actual[:, :m]), atol=3e-11)
            weights = r*(R+r*np.cos(theta))*u**6
            record(rows, interactions, label, n, h, inner, ix, weights[ix], pieces, actual,
                   wall=True)


def physical(rows, interactions, prior):
    coupled = ROOT / 'research/run43_coupled_recovery/results/fields'
    audits = {'coupled': json.loads((coupled/'independent-audit.json').read_text()),
              'hybrid': json.loads((prior/'hybrid/independent-audit.json').read_text())}
    with ZipFile(ARCHIVE) as archive:
        for index, n in enumerate(GRIDS):
            uref, wref = reference_fields(archive, prior/'levels', n, 4, {})
            hybrid = load(prior/'hybrid', n)
            assert digest((prior/'hybrid'/f'fields_{n}.npz').read_bytes()) == MANIFEST['hybrid_sha256'][str(n)]
            audit = Audit((n+1, n, n))
            inner = restrictions(audit)['inner_band']
            r, theta = audit.r[..., 0], audit.t[..., 0]
            a, _ = basis(r, theta)
            for label, (u, w) in (('coupled', load(coupled, n)), ('hybrid', hybrid),
                                  ('reference_4N', (uref, wref[:, ::8, ::4]))):
                frozen = {k: v[:, inner] for k, v in frozen_parts(audit, w, a).items()}
                samples = wref if label == 'reference_4N' else quartic_samples(w, 1)
                cont = continuum_parts(samples, n, inner, r[inner], theta)
                source = 8*np.pi*G*u[inner]**-4*audit.S[..., 0][:, inner]
                pieces, actual = decomposition(frozen, cont, u[inner], source)
                weights = (r*audit.q[..., 0]*u**6)[inner]
                if label in audits:
                    expected = audits[label]['rows'][index]['audit']['inner_boundary_band']
                    np.testing.assert_allclose(rms(weights, actual),
                        expected['independent_tensor_M_order2']['volume_weighted_rms'],
                        atol=1e-15, rtol=0)
                record(rows, interactions, label, n, audit.dr, inner, inner, weights, pieces, actual)


def deficits(rows, interactions):
    by_level = {(r['field'], r['n'], r['region']): r for r in rows}
    budgets = {}
    for item in interactions:
        if item['region'] in ('first', 'remaining'):
            budgets.setdefault((item['field'], item['n'], item['region']), []).append(item)
    for label in ('hybrid', 'coupled', 'reference_4N'):
        for coarse, fine in zip(GRIDS, GRIDS[1:]):
            def energy(n, region):
                return by_level[label, n, region]['frozen_energy_band_normalized']

            total = energy(fine, 'inner')-energy(coarse, 'inner')/8
            order = .5*np.log2(energy(coarse, 'inner')/energy(fine, 'inner'))
            for region in ('first', 'remaining'):
                change = energy(fine, region)-energy(coarse, region)/8
                low, high = budgets[label, coarse, region], budgets[label, fine, region]
                deltas = {}
                for a, b in zip(low, high):
                    assert (a['a'], a['b']) == (b['a'], b['b'])
                    delta = b['energy_band_normalized']-a['energy_band_normalized']/8
                    for key in {a['a'], a['b']}:
                        deltas[key] = deltas.get(key, 0.) + delta/(1 if a['a'] == a['b'] else 2)
                np.testing.assert_allclose(sum(deltas.values()), change, atol=1e-18)
                yield dict(field=label, coarse=coarse, fine=fine, region=region, order=order,
                           total_deficit=total, region_deficit=change, region_share=change/total,
                           **{key+'_deficit': value for key, value in deltas.items()})


def save(out, prefix, rows, interactions):
    previous = {}
    for row in rows:
        key = row['field'], row['region']
        for quantity in ('frozen', 'rr', 'metric', 'gap', 'other_principal'):
            value = row[quantity+'_rms']
            old = previous.get(key, {}).get(quantity+'_rms', 0)
            row[quantity+'_order'] = np.log2(old/value) if old > 0 and value > 0 else ''
        previous[key] = row
    write_table(out/f'{prefix}.csv', rows)
    write_table(out/f'{prefix}_interactions.csv', interactions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prior', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows, interactions = [], []
    manufactured(rows, interactions)
    save(args.output, 'manufactured', rows, interactions)
    rows, interactions = [], []
    physical(rows, interactions, args.prior)
    save(args.output, 'physical', rows, interactions)
    write_table(args.output/'deficits.csv', list(deficits(rows, interactions)))


if __name__ == '__main__':
    main()
