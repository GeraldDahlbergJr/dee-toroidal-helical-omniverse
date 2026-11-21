"""Measure the frozen boundary defect without changing any audited fields."""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from zipfile import ZipFile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for folder in ("run43_coupled_recovery", "run43_frozen_momentum"):
    sys.path.insert(0, str(ROOT / "research" / folder))
from diagnose import (ARCHIVE, GRIDS, MANIFEST, PRIOR, Audit, G, R, RMIN, RMAX,
                      basis, continuum_parts, contract, digest, dot, frozen_parts,
                      longitude, radial_matrix, reference_fields, restrictions,
                      rms, tensor, write_table)
from evaluate import load, quartic_samples
from operator_eval import divergence

BASE = "a86b8f046908f852722d1d034f4e0711d5fc9466"


def regions(band):
    return {"first_row": band[:1], "remaining_band": band[1:], "full_band": band}


def manufactured():
    rows = []
    length = RMAX-RMIN
    for n in (*GRIDS, 128, 256, 512, 1024):
        h, dt = length/n, 2*np.pi/n
        r = np.linspace(RMIN, RMAX, n+1)[:, None]
        t = (dt*np.arange(n))[None, :]
        a, _ = basis(r, t)
        # Same operators as Audit, without allocating the full 3-D geometry.
        op = SimpleNamespace(drad={2: radial_matrix(n+1, h)}, dp=dt, dt=dt)
        op.angular_derivative = Audit.angular_derivative.__get__(op)
        weight = r*(R+r*np.cos(t))*(1+.05*np.cos(t))**6
        band = np.flatnonzero((r[:, 0] > RMIN) & (r[:, 0] <= .10+1e-12))
        g = np.concatenate((np.cos(2*t), np.sin(t), 1+.2*np.cos(t)))[:, None]
        gt = np.concatenate((-2*np.sin(2*t), np.cos(t), -.2*np.sin(t)))[:, None]
        gtt = np.concatenate((-4*np.cos(2*t), -np.sin(t), -.2*np.cos(t)))[:, None]
        x = (r-RMIN)/length
        scale = (1+.05*np.cos(t))**-10
        for case in ("sine", "quadratic"):
            if case == "sine":
                k = np.pi/length
                f = np.sin(np.pi*x)
                f[[0, -1]] = 0
                fp, fpp = k*np.cos(np.pi*x), -k*k*f
                fppp = -k*k*fp
            else:
                f, fp, fpp = x*(1-x), (1-2*x)/length, -2/length**2
                fppp = np.zeros_like(x)
            w = f*g
            parts = frozen_parts(op, w, a)
            exact = divergence(w, fp*g, fpp*g, f*gt, fp*gt, f*gtt, r, t, R)
            error = scale*(sum(parts.values())-exact)
            rr = scale*(parts["rr"]-contract(tensor(fpp*g, a[0]), a[0]))
            leading = scale*h/4*contract(tensor(fppp*g, a[0]), a[0])
            if n <= 32:
                full = Audit((n+1, n, n))
                np.testing.assert_allclose(longitude(sum(parts.values()), full.p[0, 0]),
                                           full.momentum(np.ones(full.shape),
                                                         longitude(w, full.p[0, 0]), 2)
                                           +8*np.pi*G*full.S, atol=1e-8, rtol=1e-13)
            for region, ix in {"wall": np.array([0]), **regions(band)}.items():
                predicted = leading[:, ix]*((ix == 1) + 3*(ix == 0))[:, None]
                local = weight[ix]
                rows.append(dict(n=n, case=case, region=region,
                                 rms=rms(local, error[:, ix]),
                                 rr_rms=rms(local, rr[:, ix]),
                                 non_rr_rms=rms(local, (error-rr)[:, ix]),
                                 leading_rms=rms(local, predicted),
                                 remainder_rms=rms(local, error[:, ix]-predicted),
                                 first_row_weight_fraction=float(local[ix == 1].sum()/local.sum())))
    return rows


def verified_inputs(work):
    hashes = {}

    def verify(path, expected):
        actual = digest(path.read_bytes())
        assert actual == expected, path
        hashes[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT)
               else str(path.relative_to(work))] = actual

    for name, expected in MANIFEST["archive_and_code_sha256"].items():
        verify(ROOT/name, expected)
    for folder, key in ((PRIOR.parent, "new_code_sha256"), (PRIOR, "output_sha256")):
        for name, expected in MANIFEST[key].items():
            verify(folder/name, expected)
    checkpoint = ROOT/"reproduce/frozen/run42"
    manifest = json.loads((checkpoint/"manifest.json").read_text())
    artifact = checkpoint/"artifacts.zip"
    verify(artifact, manifest["source"]["artifact_digest"].removeprefix("sha256:"))
    with ZipFile(artifact) as archive:
        assert archive.testzip() is None
        for name, expected in manifest["archived_files_sha256"].items():
            actual = digest(archive.read(name))
            assert actual == expected, name
            hashes[f"{artifact.relative_to(ROOT)}:{name}"] = actual
    hybrid = work/"prior/hybrid/independent-audit.json"
    verify(hybrid, MANIFEST["full_audit_sha256"])
    archive_audit = ROOT/"research/run43_coupled_recovery/results/fields/independent-audit.json"
    verify(work/"coupled/independent-audit.json", digest(archive_audit.read_bytes()))
    for n in GRIDS:
        name = f"fields_{n}.npz"
        verify(work/"prior/hybrid"/name, MANIFEST["hybrid_sha256"][str(n)])
        expected = json.loads(archive_audit.with_name(f"result_{n}.json").read_text())["sha256"]
        verify(work/"coupled"/name, expected)
    return hashes


def archived(work):
    rows, budgets = [], []
    histories = json.loads((ARCHIVE.parent/"histories.json").read_text())
    hashes = {f'n{v["n"]}-r{v["nr"]}': v["sha256"] for v in histories}
    with ZipFile(ARCHIVE) as archive:
        assert archive.testzip() is None
        for key, expected in hashes.items():
            assert digest(archive.read(key+"/reduced.npz")) == expected, key
        for n in GRIDS:
            audit = Audit((n+1, n, n))
            r, t = audit.r[..., 0], audit.t[..., 0]
            a, _ = basis(r, t)
            band = restrictions(audit)["inner_band"]
            fields = []
            for name, folder in (("hybrid", work/"prior/hybrid"), ("coupled", work/"coupled")):
                u, w = load(folder, n)
                fields.append((name, u, w, quartic_samples(w, 1)))
            for factor in (1, 4):
                u, samples = reference_fields(archive, work/"prior/levels", n, factor, hashes)
                fields.append((f"reference_{factor}N", u, samples[:, ::8, ::factor], samples))
            for name, u, w, samples in fields:
                weight = (r*audit.q[..., 0]*u**6)[band]
                scale = u[band]**-10
                source = 8*np.pi*G*u[band]**-4*audit.S[..., 0][:, band]
                f = {k: scale*v[:, band] for k, v in frozen_parts(audit, w, a).items()}
                c = {k: scale*v for k, v in continuum_parts(samples, n, band, r[band], t).items()}
                residual = sum(f.values())-source
                continuum = sum(c.values())-source
                full_u = np.broadcast_to(u[..., None], audit.shape)
                full_w = longitude(w, audit.p[0, 0])
                original = audit.momentum(full_u, full_w, 2)
                np.testing.assert_allclose(longitude(residual, audit.p[0, 0]),
                                           original[:, band], atol=2e-14, rtol=0)
                gap = residual-continuum
                rr = f["rr"]-c["rr"]
                first = band == 1
                mixed = f["rt_metric"]-c["rt_metric"]+f["tr_metric"]-c["tr_metric"]
                tt = f["tt_metric"]-c["tt_metric"]
                rr_first, mixed_first, tt_first = rr*first[:, None], mixed*first[:, None], tt*first[:, None]
                for region, ix in regions(np.arange(len(band))).items():
                    local = weight[ix]
                    rows.append(dict(n=n, field=name, region=region,
                                     rms=rms(local, residual[:, ix]),
                                     continuum_rms=rms(local, continuum[:, ix]),
                                     gap_rms=rms(local, gap[:, ix]),
                                     rr_gap_rms=rms(local, rr[:, ix]),
                                     non_rr_gap_rms=rms(local, (gap-rr)[:, ix]),
                                     rt_tr_metric_gap_rms=rms(local, mixed[:, ix]),
                                     tt_metric_gap_rms=rms(local, tt[:, ix]),
                                     band_weight_fraction=float(local.sum()/weight.sum())))
                terms = {"rr_first": rr_first, "rr_remaining": rr-rr_first,
                         "rt_tr_metric_first": mixed_first, "rt_tr_metric_remaining": mixed-mixed_first,
                         "tt_metric_first": tt_first, "tt_metric_remaining": tt-tt_first,
                         "other_first": gap*first[:, None]-rr_first-mixed_first-tt_first,
                         "other_remaining": (gap-rr-mixed-tt)*(~first)[:, None],
                         "continuum": continuum}
                np.testing.assert_allclose(sum(terms.values()), residual, atol=1e-16)
                energy = dot(weight, residual, residual)
                np.testing.assert_allclose(sum(dot(weight, v, residual) for v in terms.values()),
                                           energy, atol=1e-19)
                terms["all_first_gap"] = gap*first[:, None]
                for term, defect in terms.items():
                    budgets.append(dict(n=n, field=name, term=term, energy=energy,
                                        allocation=dot(weight, defect, residual),
                                        removed_energy=dot(weight, residual-defect, residual-defect)))
    return rows, budgets


def deficits(budgets):
    previous = {}
    for row in budgets:
        key = row["field"], row["term"]
        if key in previous:
            coarse = previous[key]
            deficit = row["energy"]-coarse["energy"]/8
            attributed = row["allocation"]-coarse["allocation"]/8
            removed = row["removed_energy"]-coarse["removed_energy"]/8
            yield {**row, "coarse": coarse["n"], "deficit": deficit, "signed_share": attributed/deficit,
                   "allocated_deficit": attributed, "deficit_after_removal": removed,
                   "order_after_removal": float(np.log2(coarse["removed_energy"]/row["removed_energy"])/2)}
        previous[key] = row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    hashes = verified_inputs(args.work)
    physical, budgets = archived(args.work)
    analytic = manufactured()
    args.out.mkdir(exist_ok=False)
    write_table(args.out/"manufactured.csv", analytic, ("case", "region"))
    write_table(args.out/"archived.csv", physical, ("field", "region"))
    write_table(args.out/"deficits.csv", list(deficits(budgets)))
    for name, folder in (("hybrid", "prior/hybrid"), ("coupled", "coupled")):
        report = json.loads((args.work/folder/"independent-audit.json").read_text())
        for record in report["rows"]:
            measured = next(v["rms"] for v in physical if v["n"] == record["n"]
                            and v["field"] == name and v["region"] == "full_band")
            expected = record["audit"]["inner_boundary_band"]["independent_tensor_M_order2"]["volume_weighted_rms"]
            np.testing.assert_allclose(measured, expected, atol=2e-14, rtol=0)
        del report["rows"]
        (args.out/f"{name}-acceptance.json").write_text(json.dumps(report, indent=2)+"\n")
    (args.out/"provenance.json").write_text(json.dumps(
        {"base_commit": BASE, "python": sys.version.split()[0], "numpy": np.__version__,
         "verified_sha256": hashes,
         "reference_members_sha256": {f'n{v["n"]}-r{v["nr"]}': v["sha256"] for v in
                                      json.loads((ARCHIVE.parent/"histories.json").read_text())},
         "reference_4N_sha256": {str(n): digest(
             (args.work/"prior/levels"/f"a{4*n}-r{8*n}"/"reduced.npz").read_bytes()) for n in GRIDS}},
        indent=2)+"\n")


if __name__ == "__main__":
    main()
