"""Decompose the Run 43 tensor operator and test radial substitutions."""
import argparse
import csv
import hashlib
from io import BytesIO
from itertools import product
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
from numpy.polynomial import Polynomial

ROOT = Path(__file__).resolve().parents[2]
for name in ("reproduce", "research/run43_asymptotic"):
    sys.path.insert(0, str(ROOT / name))
from diagnose_toroidal_boundary_convergence import Audit, R, RMIN, RMAX, G, radial_matrix
from operator_eval import J, longitude, traces

AXES = "rtp"
PAIRS = tuple(product(range(3), repeat=2))
GRIDS = (8, 16, 32, 64)
ARCHIVE = ROOT / "research/run43_reference/results/Run43_Radial_Reference_Fields.zip"
PRIOR = ROOT / "research/run43_asymptotic/results"
MANIFEST = json.loads((PRIOR / "provenance.json").read_text())


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def radial(matrix, value):
    return np.einsum("ij,...jt->...it", matrix, value)


def second(n, h, count):
    result = np.zeros((n, n))
    for row in range(n):
        start = min(max(row-count//2, 0), n-count)
        cols = np.arange(start, start+count)
        offsets = (cols-row).astype(float)
        rhs = np.zeros(count)
        rhs[2] = 2
        result[row, cols] = np.linalg.solve(offsets**np.arange(count)[:, None], rhs)/h**2
    return result


def multiply(matrix, value):
    return np.einsum("ij,j...->i...", matrix, value)


def tensor(v, a):
    gradient = v[:, None]*a[None, :]
    trace = np.einsum("ii...->...", gradient)
    return gradient + gradient.swapaxes(0, 1) - (2/3)*np.eye(3)[:, :, None, None]*trace


def contract(t, b):
    return np.einsum("ijrt,jrt->irt", t, b)


def basis(r, theta):
    c, s = np.cos(theta), np.sin(theta)
    er = np.stack(np.broadcast_arrays(c, 0*c, s))
    et = np.stack(np.broadcast_arrays(-s, 0*c, c))
    ep = np.stack(np.broadcast_arrays(0*c, 1+0*c, 0*c))
    q = R+r*c
    a = (er+0*r, et/r, ep/q)
    zero = np.zeros_like(a[0])
    derivatives = ((zero, -et/r**2, -c*ep/q**2),
                   (et+0*r, -er/r, r*s*ep/q**2),
                   tuple(multiply(J, v) for v in a))
    return a, derivatives


def frozen_parts(audit, w, a):
    # Reduce at phi=0 by equivariance. A Cartesian vector rotates by Q,
    # and a Cartesian tensor rotates by Q T Q^T.
    dr = audit.drad[2]
    dp = np.sin(audit.dp)/audit.dp
    q = np.array([[np.cos(audit.dp), -np.sin(audit.dp), 0],
                  [np.sin(audit.dp), np.cos(audit.dp), 0], [0, 0, 1]])

    def dt(value):
        return audit.angular_derivative(value, -1, audit.dt, 2)

    def phi_tensor(t):
        return (np.einsum("ia,abrt,jb->ijrt", q, t, q) -
                np.einsum("ia,abrt,jb->ijrt", q.T, t, q.T))/(2*audit.dp)

    derivatives = (lambda v: radial(dr, v), dt, lambda v: dp*multiply(J, v))
    first = [derive(w) for derive in derivatives]
    tensor_derivatives = (*derivatives[:2], phi_tensor)
    parts = {}
    for i, j in PAIRS:
        name = AXES[i]+AXES[j]
        principal = contract(tensor(derivatives[i](first[j]), a[j]), a[i])
        complete = contract(tensor_derivatives[i](tensor(first[j], a[j])), a[i])
        parts[name] = principal
        parts[name+"_metric"] = complete-principal
    return parts


def continuum_parts(w, n, indices, r, theta):
    """Average the quartic left/right jets and differentiate each product."""
    radial_intervals, angles = w.shape[1]-1, w.shape[2]
    angular = slice(None, None, angles//n)
    nodes = indices*(radial_intervals//n)
    k = np.fft.fftfreq(angles, d=1/angles)
    first = 1j*k
    first[angles//2] = 0
    spectrum = np.fft.fft(w, axis=-1)
    wt = np.fft.ifft(spectrum*first, axis=-1).real
    wtt = np.fft.ifft(-spectrum*k**2, axis=-1).real
    grid = np.linspace(0, 1, 5)
    polys = [Polynomial.fromroots(np.delete(grid, j)) for j in range(5)]
    polys = [p/p(x) for p, x in zip(polys, grid)]
    scale = radial_intervals/(4*(RMAX-RMIN))

    def trace_mean(field, order):
        values = []
        for endpoint in (1, 0):
            element = nodes[:, None]+np.arange(5)-4*endpoint
            coeff = scale**order*np.array([p.deriv(order)(endpoint) for p in polys])
            samples = field[:, element, angular]-field[:, nodes, None, angular]
            values.append(np.einsum("j,irjt->irt", coeff, samples))
        return (values[0]+values[1])/2

    v = w[:, nodes, angular]
    vr, vt = trace_mean(w, 1), wt[:, nodes, angular]
    vp = multiply(J, v)
    vrt, vrp, vtp = trace_mean(wt, 1), multiply(J, vr), multiply(J, vt)
    first = (vr, vt, vp)
    second = ((trace_mean(w, 2), vrt, vrp),
              (vrt, wtt[:, nodes, angular], vtp),
              (vrp, vtp, multiply(J, vp)))
    a, da = basis(r, theta)
    parts = {}
    for i, j in PAIRS:
        name = AXES[i]+AXES[j]
        parts[name] = contract(tensor(second[i][j], a[j]), a[i])
        parts[name+"_metric"] = contract(tensor(first[j], da[i][j]), a[i])
    expected = sum(traces(w, n, indices, r, theta, R, RMAX-RMIN))/2
    np.testing.assert_allclose(sum(parts.values()), expected, rtol=0, atol=2e-14)
    return parts


def corrections(audit, w, a, da, frozen):
    d1 = audit.drad[2]
    composed = d1@d1
    compact = second(len(d1), audit.dr, 3)
    fourth_wall = second(len(d1), audit.dr, 6)
    deltas = {"compact_rr": compact-composed,
              "wall_rr_compact": compact-composed,
              "wall_rr_fourth": fourth_wall-composed,
              "fourth_rr": fourth_wall-composed}
    for name in ("wall_rr_compact", "wall_rr_fourth"):
        deltas[name][2:-2] = 0
    fixes = {name: contract(tensor(radial(delta, w), a[0]), a[0])
             for name, delta in deltas.items()}
    # Product rule diagnostic using the same frozen nodal first derivatives.
    first = (radial(d1, w), audit.angular_derivative(w, -1, audit.dt, 2),
             np.sin(audit.dp)/audit.dp*multiply(J, w))
    metric = sum(contract(tensor(first[j], da[i][j]), a[i]) -
                 frozen[AXES[i]+AXES[j]+"_metric"] for i, j in PAIRS)
    return {**fixes, "product_rule_metric": metric,
            "compact_rr_and_metric": fixes["compact_rr"]+metric}


def rms(weight, value):
    return float(np.sqrt(np.sum(weight*np.sum(value**2, axis=0))/weight.sum()))


def write_table(path, rows, series=None, level="n", value="rms"):
    if series:
        prior = {}
        for row in rows:
            key = tuple(row[k] for k in series)
            previous = prior.get(key)
            row["order"] = (np.log2(previous[value]/row[value])
                            if previous and previous[value] and row[value] else "")
            assert not previous or previous[level]*2 == row[level]
            prior[key] = row
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def stencil_tables(out):
    coefficients, errors = [], []
    for n in (8, 16, 32, 64, 128, 256):
        h = 1/n
        d1 = radial_matrix(n+1, h)
        methods = {"composed": d1@d1, "compact": second(n+1, h, 3),
                   "fourth": second(n+1, h, 6)}
        x = np.linspace(0, 1, n+1)
        for name, matrix in methods.items():
            for row in range(5):
                coeff = h**2*matrix[row]
                offsets = np.arange(n+1)-row
                if n == 8:
                    coefficients.append(dict(method=name, row=row, **{
                        f"offset_{k}": coeff[row+k] if 0 <= row+k <= n else 0
                        for k in range(-3, 6)},
                        cubic_moment=coeff@offsets**3,
                        quartic_moment=coeff@offsets**4))
                np.testing.assert_allclose([coeff@offsets**k for k in range(3)],
                                           [0, 0, 2], atol=1e-14)
                for label, f, exact in (("cubic", x**3, 6*x),
                                         ("exponential", np.exp(x), np.exp(x)),
                                         ("sine", np.sin(2*x), -4*np.sin(2*x))):
                    errors.append(dict(n=n, method=name, row=row, field=label,
                                       error=abs(matrix[row]@(f-f[row])-exact[row])))
    write_table(out/"stencils.csv", coefficients)
    write_table(out/"manufactured.csv", errors, ("method", "row", "field"), value="error")


def reference_fields(archive, levels, n, factor, hashes):
    radial_count, angles = 8*n, factor*n
    if factor == 1:
        key = f"n{n}-r{radial_count}"
        raw = archive.read(key+"/reduced.npz")
        expected = hashes[key]
    else:
        raw = (levels/f"a{angles}-r{radial_count}"/"reduced.npz").read_bytes()
        expected = REFINED[angles, radial_count]
    assert digest(raw) == expected
    with np.load(BytesIO(raw)) as data:
        u, w = data["psi"], data["W"]
    assert u.shape == (radial_count+1, angles) and w.shape == (3, *u.shape)
    assert np.all(u[[0, -1]] == 1) and np.count_nonzero(w[:, [0, -1]]) == 0
    return u[::8, ::factor], w


REFINED = {(row["angles"], row["radial"]): row["sha256"]
           for row in json.loads((PRIOR/"histories.json").read_text()) if "angles" in row}


def restrictions(audit):
    r = audit.r[:, 0, 0]
    interior = (r > r[0]) & (r < r[-1])
    return {"first_row": np.array([1]), "inner_band": np.flatnonzero(interior & (r <= .10+1e-12)),
            "outer_band": np.flatnonzero(interior & (r >= .20-1e-12))}


def axis_groups(parts):
    return {**{f"outer_{axis}": sum(v for k, v in parts.items() if k[0] == axis)
               for axis in AXES},
            "rt_tr_metric": parts["rt_metric"]+parts["tr_metric"],
            "poloidal_metric": parts["rt_metric"]+parts["tr_metric"]+parts["tt_metric"],
            "mixed_principal": sum(v for k, v in parts.items() if len(k) == 2 and k[0] != k[1]),
            "mixed": sum(v for k, v in parts.items() if k[0] != k[1]),
            "metric": sum(v for k, v in parts.items() if k.endswith("_metric")),
            "all_principal": sum(v for k, v in parts.items() if len(k) == 2)}


def dot(weight, x, y):
    return float(np.sum(weight*np.sum(x*y, axis=0))/weight.sum())


def diagnose_reference(archive, hashes, levels, n, factor, audit, weight,
                       hybrid_residual, terms, comparisons, allocations):
    u, ref = reference_fields(archive, levels, n, factor, hashes)
    w = ref[:, ::8, ::factor]
    r, theta = audit.r[:, :, 0], audit.t[:, :, 0]
    a, da = basis(r, theta)
    frozen = frozen_parts(audit, w, a)
    phi = audit.p[0, 0]
    full_u = np.broadcast_to(u[..., None], audit.shape)
    source = 8*np.pi*G*u**-4*audit.S[..., 0]
    residual = u**-10*sum(frozen.values())-source
    np.testing.assert_allclose(longitude(residual, phi),
                               audit.momentum(full_u, longitude(w, phi), 2),
                               atol=2e-14, rtol=0)
    fixes = corrections(audit, w, a, da, frozen)
    z = u**-10
    for region, ix in restrictions(audit).items():
        f = {name: z[ix]*v[:, ix] for name, v in frozen.items()}
        c = {name: z[ix]*v for name, v in continuum_parts(ref, n, ix, r[ix], theta).items()}
        f.update(axis_groups(f))
        c.update(axis_groups(c))
        gaps = {name: f[name]-c[name] for name in f}
        # Resolve the tensor's rr principal contribution into its Cartesian
        # vector Laplacian and one-third gradient-of-divergence pieces.
        normal = np.sum(a[0][:, ix]*gaps["rr"], axis=0)*a[0][:, ix]
        gaps["rr_grad_div"] = normal/4
        gaps["rr_laplacian"] = gaps["rr"]-normal/4
        frozen_residual = residual[:, ix]
        continuum_residual = c["all_principal"]+c["metric"]-source[:, ix]
        gap = frozen_residual-continuum_residual
        np.testing.assert_allclose(gap, gaps["all_principal"]+gaps["metric"], atol=1e-15)
        g = weight[ix]
        denom = dot(g, gap, gap)
        for name, delta in gaps.items():
            terms.append(dict(n=n, factor=factor, region=region, term=name,
                              frozen_rms=rms(g, f[name]) if name in f else "",
                              continuum_rms=rms(g, c[name]) if name in c else "",
                              rms=rms(g, delta), signed_projection=dot(g, delta, gap),
                              gap_projection_share=dot(g, delta, gap)/denom,
                              gap_after_substitution=rms(g, gap-delta),
                              residual_after_substitution=rms(g, frozen_residual-delta)))
        h = hybrid_residual[..., 0][:, ix]
        rr_wall = gaps["rr"]*(ix == 1)[None, :, None]
        attribution = {
            "rr_first_row": rr_wall,
            "rr_other_rows": gaps["rr"]-rr_wall,
            "rt_tr_metric": gaps["rt_tr_metric"],
            "tt_metric": gaps["tt_metric"],
            "tt_principal": gaps["tt"],
            "other_operator": gap-gaps["rr"]-gaps["poloidal_metric"]-gaps["tt"],
            "continuum_residual": continuum_residual,
            "hybrid_minus_reference_frozen": h-frozen_residual,
        }
        np.testing.assert_allclose(sum(attribution.values()), h, atol=2e-16, rtol=0)
        assert abs(sum(dot(g, v, h) for v in attribution.values())-dot(g, h, h)) < 1e-18
        for name, v in attribution.items():
            allocations.append(dict(n=n, factor=factor, region=region, term=name,
                                    hybrid_energy=dot(g, h, h), allocation=dot(g, v, h)))
        candidates = {"frozen": frozen_residual, "continuum": continuum_residual,
                      **{name: frozen_residual+z[ix]*v[:, ix] for name, v in fixes.items()}}
        for name, v in candidates.items():
            comparisons.append(dict(n=n, factor=factor, region=region, candidate=name,
                                    rms=rms(g, v), gap_rms=rms(g, v-continuum_residual)))


def diagnose_hybrid(audit, u, w, summary, rows):
    a, da = basis(audit.r[:, :, 0], audit.t[:, :, 0])
    reduced = w[..., 0]
    fixes = corrections(audit, reduced, a, da, frozen_parts(audit, reduced, a))
    original = audit.momentum(u, w, 2)
    phi = audit.p[0, 0]
    for band, ix in restrictions(audit).items():
        if band == "first_row":
            continue
        key = band.replace("_band", "_boundary_band")
        mask = np.zeros(audit.shape, dtype=bool)
        mask[ix] = True
        for name, delta in {"frozen": np.zeros_like(w[..., 0]), **fixes}.items():
            corrected = original+u**-10*longitude(delta, phi)
            values = {"momentum": audit.norms(corrected, u, mask)["volume_weighted_rms"],
                      "hamiltonian": summary[key]["independent_flux_H"]["volume_weighted_rms"]}
            for metric, value in values.items():
                rows.append(dict(n=audit.shape[1], band=band, candidate=name, metric=metric,
                                 rms=value))
    return original


def acceptance(rows, frozen_report):
    checks = []
    for row in rows:
        if row["order"] == "":
            continue
        check = dict(candidate=row["candidate"], band=row["band"], metric=row["metric"],
                     pair=[row["n"]//2, row["n"]], order=row["order"],
                     passed=bool(np.isfinite(row["order"]) and row["order"] >= 1.5))
        checks.append(check)
    frozen = [check for check in checks if check["candidate"] == "frozen"]
    assert len(frozen) == len(frozen_report["checks"]) == 12
    metric = {"independent_flux_H": "hamiltonian", "independent_tensor_M_order2": "momentum"}
    index = {(c["band"], c["metric"], tuple(c["pair"])): c for c in frozen}
    for original in frozen_report["checks"]:
        key = (original["band"].replace("_boundary", ""), metric[original["metric"]],
               tuple(original["pair"]))
        assert original["passed"] == index[key]["passed"]
        np.testing.assert_allclose(original["order"], index[key]["order"], atol=1e-14, rtol=0)
    scores = {}
    for candidate in dict.fromkeys(c["candidate"] for c in checks):
        subset = [c for c in checks if c["candidate"] == candidate]
        assert len(subset) == 12
        scores[candidate] = dict(passed=sum(c["passed"] for c in subset), total=12)
    return dict(scores=scores, checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hybrid", type=Path)
    parser.add_argument("levels", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--factors", type=int, nargs="+", default=[1, 4])
    args = parser.parse_args()
    for name, expected in MANIFEST["archive_and_code_sha256"].items():
        assert digest((ROOT/name).read_bytes()) == expected, name
    assert digest((ROOT/"research/run43_asymptotic/operator_eval.py").read_bytes()) == (
        MANIFEST["new_code_sha256"]["operator_eval.py"])
    full_audit = (args.hybrid/"independent-audit.json").read_bytes()
    assert digest(full_audit) == MANIFEST["full_audit_sha256"]
    report = json.loads(full_audit)
    terms, comparisons, hybrid, allocations = [], [], [], []
    with ZipFile(ARCHIVE) as archive:
        hashes = {f'n{row["n"]}-r{row["nr"]}': row["sha256"] for row in
                  json.loads((ROOT/"research/run43_reference/results/histories.json").read_text())}
        for n, record in zip(GRIDS, report["rows"]):
            raw = (args.hybrid/f"fields_{n}.npz").read_bytes()
            assert digest(raw) == MANIFEST["hybrid_sha256"][str(n)] == record["sha256"]
            with np.load(BytesIO(raw)) as data:
                u, w = data["psi"], data["W_cartesian"]
            audit = Audit(u.shape)
            hybrid_residual = diagnose_hybrid(audit, u, w, record["audit"], hybrid)
            weight = (audit.r*audit.q*u**6)[..., 0]
            for factor in args.factors:
                diagnose_reference(archive, hashes, args.levels, n, factor, audit, weight,
                                   hybrid_residual, terms, comparisons, allocations)
    with (PRIOR/"residuals.csv").open() as handle:
        prior = {(int(row["n"]), int(row["factor"]), row["region"]): row
                 for row in csv.DictReader(handle)}
    for row in comparisons:
        key = (row["n"], row["factor"], row["region"])
        if key not in prior:
            continue
        old = prior[key]
        if row["candidate"] in ("frozen", "continuum"):
            field = {"frozen": "frozen", "continuum": "independent"}[row["candidate"]]
            np.testing.assert_allclose(row["rms"], float(old[field]), atol=2e-14, rtol=0)
        if row["candidate"] == "frozen":
            np.testing.assert_allclose(row["gap_rms"], float(old["operator_gap"]),
                                       atol=2e-14, rtol=0)
    args.out.mkdir(exist_ok=False)
    stencil_tables(args.out)
    write_table(args.out/"terms.csv", terms, ("factor", "region", "term"))
    write_table(args.out/"reference.csv", comparisons, ("factor", "region", "candidate"))
    write_table(args.out/"hybrid.csv", hybrid, ("band", "candidate", "metric"))
    previous, deficits = {}, []
    for row in allocations:
        key = (row["factor"], row["region"], row["term"])
        if key in previous:
            coarse = previous[key]
            deficit = row["hybrid_energy"]-coarse["hybrid_energy"]/8
            allocation = row["allocation"]-coarse["allocation"]/8
            deficits.append({**row, "coarse": coarse["n"], "hybrid_deficit": deficit,
                             "allocated_deficit": allocation, "share": allocation/deficit})
        previous[key] = row
    write_table(args.out/"deficits.csv", deficits)
    (args.out/"acceptance.json").write_text(json.dumps(acceptance(hybrid, report), indent=2)+"\n")
    del report["rows"]
    (args.out/"frozen-acceptance.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
