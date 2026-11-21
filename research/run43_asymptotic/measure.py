"""Compare both operators and decompose the unchanged hybrid momentum deficit."""
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from importlib.metadata import version
from io import BytesIO
from itertools import combinations_with_replacement
from pathlib import Path
from zipfile import ZipFile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run43_angular"))
from analyze import ARCHIVE, Audit, GRIDS, ROOT, digest, restrict
from operator_eval import divergence, longitude, traces
from solve import reference

BASE = "6e0e9f867a4c0f27b34a0d930c3af336e984ff8b"
OLD = ROOT / "research/run43_angular/results"
MAX_FACTOR = {8: 64, 16: 64, 32: 32, 64: 8}


def check_operator(audit):
    r, theta = audit.r[1:-1, :, 0], audit.t[:, :, 0]
    c, s = np.cos(theta), np.sin(theta)
    major = float(audit.q[0, 0, 0]-audit.r[0, 0, 0])
    q = major+r*c
    # Cartesian W=(0,0,z²) and cylindrical W=q² e_rho both have
    # closed-form div(LW) from Laplacian(W)+grad(div(W))/3.
    for axis, x, xr, xt, xrt, xtt, expected in (
        (2, r*s, s, r*c, c, -r*s, 8/3),
        (0, q, c, -r*s, -s, -r*c, 4),
    ):
        direction = np.eye(3)[axis, :, None, None]
        actual = divergence(*(direction*v for v in
                            (x*x, 2*x*xr, 2*xr*xr, 2*x*xt,
                             2*(xt*xr+x*xrt), 2*(xt*xt+x*xtt))), r, theta, major)
        np.testing.assert_allclose(actual, np.broadcast_to(direction*expected, actual.shape),
                                   atol=2e-12, rtol=0)


def inner(weight, a, b):
    return float(np.sum(weight*np.sum(a*b, axis=0))/weight.sum())


def norm(weight, a):
    return np.sqrt(inner(weight, a, a))


def contributions(weight, parts, hybrid):
    terms = {f"{a}*{b}": (1 if a == b else 2)*inner(weight, parts[a], parts[b])
             for a, b in combinations_with_replacement(parts, 2)}
    np.testing.assert_allclose(sum(terms.values()), inner(weight, hybrid, hybrid), rtol=2e-14)
    return terms


def measure_grid(n, archive, data, hybrid, provenance, histories, factors, radial_factor=8):
    """Yield residual vectors on matched audit nodes, with common hybrid weights."""
    path = hybrid / f"fields_{n}.npz"
    assert digest(path) == provenance["verified_hybrid_sha256"][str(n)]
    with np.load(path) as fields:
        psi, w = fields["psi"], fields["W_cartesian"]
    audit = Audit(psi.shape)
    check_operator(audit)
    indices = np.flatnonzero((audit.r[:, 0, 0] <= .10+1e-12) & (np.arange(n+1) > 0))
    h = audit.momentum(psi, w, 2)[:, indices]
    weight = (audit.r*audit.q*psi**6)[indices]
    r, theta = audit.r[indices, :, 0], audit.t[:, :, 0]
    major = float(audit.q[0, 0, 0]-audit.r[0, 0, 0])
    width = float(audit.r[-1, 0, 0]-audit.r[0, 0, 0])
    radial = radial_factor*n
    for factor in factors:
        angles = factor*n
        if factor == 1:
            raw = archive.read(f"n{n}-r{radial}/reduced.npz")
            assert hashlib.sha256(raw).hexdigest() == provenance["verified_radial_sha256"][f"n{n}-r{radial}"]
            path = BytesIO(raw)
        else:
            folder = data / f"a{angles}-r{radial}"
            info = json.loads((folder / "result.json").read_text())
            assert (angles, radial) == (info["angles"], info["radial"])
            path = folder / "reduced.npz"
            assert digest(path) == info["sha256"]
            histories.append(dict(n=n, factor=factor, **info))
        fields = restrict(path, n, angles, radial)
        frozen = audit.momentum(*fields, 2)[:, indices]
        if factor == 1:
            path.seek(0)
        with np.load(path) as reduced:
            left, right = traces(reduced["W"], n, indices, r, theta, major, width)
        u = fields[0][indices]
        source = 8*np.pi*reference.G*u**-4*audit.S[:, indices]
        phi = audit.p[0, 0]
        left = u**-10*longitude(left, phi)-source
        right = u**-10*longitude(right, phi)-source
        yield factor, weight, h, frozen, (left+right)/2, (left-right)/2


def analyze(data, hybrid, out):
    out.mkdir()
    provenance = json.loads((OLD / "provenance.json").read_text())
    for path, expected in provenance["sha256"].items():
        assert digest(ROOT / path) == expected, path
    previous = json.loads((ROOT / "research/run43_reference/results/histories.json").read_text())
    histories = [dict(row, factor=1) for row in previous if row["nr"] == 8*row["n"]]
    archived = {(int(row["n"]), row["region"]): row for row in
                csv.DictReader((OLD / "decomposition.csv").open())}
    archived_pairs = {(row["region"], int(row["coarse"]), row["term"]): row for row in
                      csv.DictReader((OLD / "pairs.csv").open())}
    archived_solves = {(row["angles"], row["radial"]): row["sha256"] for row in
                       json.loads((OLD / "histories.json").read_text())}
    rows, floors, per_grid, pairs, orders, radial_rows = [], [], {}, [], [], []
    with ZipFile(ROOT / ARCHIVE) as archive:
        for n in GRIDS:
            factors = [2**power for power in range(MAX_FACTOR[n].bit_length())]
            levels = list(measure_grid(n, archive, data, hybrid, provenance, histories, factors))
            (half_level,) = measure_grid(n, archive, data, hybrid, provenance, histories,
                                         factors=(4,), radial_factor=4)
            for region, selection in (("first_row", slice(0, 1)), ("inner_band", slice(None))):
                vectors = {}
                previous_delta = {}
                for factor, weight, h, frozen, independent, jump in levels:
                    weight = weight[selection]
                    h, frozen, independent, jump = (v[:, selection] for v in
                                                     (h, frozen, independent, jump))
                    vectors[factor] = (frozen, independent)
                    row = dict(n=n, region=region, factor=factor, angles=factor*n, radial=8*n,
                               hybrid=norm(weight, h), frozen=norm(weight, frozen),
                               independent=norm(weight, independent),
                               operator_gap=norm(weight, frozen-independent),
                               trace_half_jump=norm(weight, jump),
                               independent_left=norm(weight, independent+jump),
                               independent_right=norm(weight, independent-jump))
                    if factor <= 4:
                        old_key = {1: "archived_rms", 2: "refined_2n_rms", 4: "refined_4n_rms"}[factor]
                        np.testing.assert_allclose(row["frozen"], float(archived[n, region][old_key]),
                                                   rtol=1e-10, atol=1e-15)
                    if factor//2 in vectors and factor > 1:
                        for label, current, earlier in zip(("frozen", "independent"),
                                                           vectors[factor], vectors[factor//2]):
                            delta = norm(weight, current-earlier)
                            row[label+"_delta"] = delta
                            if label in previous_delta:
                                row[label+"_delta_order"] = np.log2(previous_delta[label]/delta)
                            previous_delta[label] = delta
                            estimate = current+(current-earlier)/3
                            floors.append(dict(n=n, region=region, factor=factor, operator=label,
                                               rms=norm(weight, estimate), correction_rms=delta/3))
                    rows.append(row)
                _, _, _, half_frozen, half_independent, half_jump = half_level
                full_frozen, full_independent = vectors[4]
                half_frozen, half_independent, half_jump = (v[:, selection] for v in
                                                           (half_frozen, half_independent, half_jump))
                radial_rows.append(dict(n=n, region=region, factor=4, radial_coarse=4*n,
                                        radial_fine=8*n,
                                        frozen_delta=norm(weight, full_frozen-half_frozen),
                                        independent_delta=norm(weight, full_independent-half_independent),
                                        independent_coarse=norm(weight, half_independent),
                                        independent_fine=norm(weight, full_independent),
                                        coarse_trace_half_jump=norm(weight, half_jump)))
                first = vectors[1][0]
                finest, independent = vectors[factor]
                angular, hybrid_gap = first-finest, h-first
                schemes = dict(
                    refined=contributions(weight, dict(reference=finest, angular=angular,
                                                        hybrid_gap=hybrid_gap), h),
                    split=contributions(weight, dict(operator=finest-independent, reference=independent,
                                                     angular=angular, hybrid_gap=hybrid_gap), h))
                reference4 = vectors[4][0]
                schemes["archived_4n"] = contributions(
                    weight, dict(reference=reference4, angular=first-reference4,
                                 hybrid_gap=hybrid_gap), h)
                per_grid[n, region] = (row, schemes)

    for info in histories:
        if info["factor"] in (2, 4):
            assert info["sha256"] == archived_solves[info["angles"], info["radial"]]

    for region in ("first_row", "inner_band"):
        for coarse, fine in zip(GRIDS, GRIDS[1:]):
            a, sa = per_grid[coarse, region]
            b, sb = per_grid[fine, region]
            identity = dict(region=region, coarse=coarse, fine=fine,
                            coarse_factor=a["factor"], fine_factor=b["factor"])
            deficit = b["hybrid"]**2-a["hybrid"]**2/8
            for scheme, terms in sa.items():
                for name, value in terms.items():
                    contribution = sb[scheme][name]-value/8
                    if scheme == "archived_4n":
                        old_term = {
                            "reference*reference": "reference_sq",
                            "angular*angular": "angular_sq",
                            "hybrid_gap*hybrid_gap": "hybrid_minus_archived_sq",
                            "reference*angular": "cross_ref_angular",
                            "reference*hybrid_gap": "cross_ref_hybrid",
                            "angular*hybrid_gap": "cross_angular_hybrid",
                        }[name]
                        old = archived_pairs[region, coarse, old_term]
                        np.testing.assert_allclose(contribution, float(old["term_gate_deficit"]),
                                                   rtol=0, atol=1e-20)
                    factors = (dict(coarse_factor=4, fine_factor=4)
                               if scheme == "archived_4n" else {})
                    pairs.append(dict(**(identity | factors), scheme=scheme, term=name,
                                      hybrid_deficit=deficit, contribution=contribution,
                                      signed_share=contribution/deficit))
                np.testing.assert_allclose(
                    sum(row["contribution"] for row in pairs[-len(terms):]), deficit,
                    rtol=0, atol=1e-20)
            for factor in sorted({r["factor"] for r in rows if r["n"] == coarse} &
                                 {r["factor"] for r in rows if r["n"] == fine}):
                x, y = (next(r for r in rows if (r["n"], r["region"], r["factor"]) == (n, region, factor))
                        for n in (coarse, fine))
                orders.append(dict(region=region, coarse=coarse, fine=fine, factor=factor,
                                   **{key: np.log2(x[key]/y[key])
                                      for key in ("hybrid", "frozen", "independent", "operator_gap")}))
    for name, table in (("residuals", rows), ("floors", floors), ("orders", orders),
                        ("pairs", pairs), ("radial_control", radial_rows)):
        if not table:
            continue
        with (out / f"{name}.csv").open("w", newline="") as stream:
            fields = list(dict.fromkeys(key for row in table for key in row))
            writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(table)
    (out / "histories.json").write_text(json.dumps(histories, indent=2)+"\n")

    run = subprocess.run([sys.executable, ROOT / "research/audit_run43_arrays.py", hybrid],
                         cwd=ROOT, capture_output=True, text=True)
    assert run.returncode == 1 and not run.stderr, run.stderr
    report_path = hybrid / "independent-audit.json"
    report = json.loads(report_path.read_text())
    assert (report["passed"], report["total"], len(report["checks"])) == (9, 12, 12)
    del report["rows"]
    compact = (json.dumps(report, indent=2)+"\n").encode()
    assert compact == (OLD / "hybrid-audit.json").read_bytes()
    (out / "hybrid-audit.json").write_bytes(compact)
    manifest = dict(base_commit=BASE, python=sys.version.split()[0], numpy=np.__version__,
                    scipy=version("scipy"), pymetis=version("pymetis"), mpmath=version("mpmath"),
                    new_code_sha256={p.name: digest(p) for p in sorted(Path(__file__).parent.glob("*.py"))}, archive_and_code_sha256=provenance["sha256"],
                    hybrid_sha256=provenance["verified_hybrid_sha256"],
                    full_audit_sha256=digest(report_path),
                    output_sha256={p.name: digest(p) for p in sorted(out.iterdir())})
    (out / "provenance.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("data", "hybrid", "out"):
        parser.add_argument(name, type=lambda p: Path(p).resolve())
    args = parser.parse_args()
    analyze(args.data, args.hybrid, args.out)
