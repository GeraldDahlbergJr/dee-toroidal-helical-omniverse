"""Compare v1, frozen and independent continuum residuals on verified arrays."""
import argparse
from importlib.metadata import version
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np

from auditor import Audit, G, R, ROOT, RMIN, RMAX, VERSION, masks, momentum
for directory in ("run43_coupled_recovery", "run43_frozen_momentum", "run43_frozen_boundary"):
    sys.path.insert(0, str(ROOT / "research" / directory))
from diagnose import ARCHIVE, GRIDS, digest, reference_fields, write_table
from evaluate import load, quartic_samples
from investigate import verified_inputs
from operator_eval import longitude, traces

BASE = "3715241456898277e290fe135222582d538c182e"


def continuum(audit, psi, samples, n, rows):
    """Use prior independent tensor traces, differentiating before restriction."""
    r, theta = audit.r[rows, :, 0], audit.t[..., 0]
    divergence = sum(traces(samples, n, rows, r, theta, R, RMAX-RMIN))/2
    return psi[rows]**-10*longitude(divergence, audit.p[0, 0]) - (
        8*np.pi*G*psi[rows]**-4*audit.S[:, rows])


def measure(work):
    results = []
    histories = json.loads((ARCHIVE.parent/"histories.json").read_text())
    hashes = {f'n{v["n"]}-r{v["nr"]}': v["sha256"] for v in histories}
    with ZipFile(ARCHIVE) as archive:
        for n in GRIDS:
            audit = Audit((n+1, n, n))
            azimuth = audit.p[0, 0]
            fields = []
            for name, folder in (("hybrid", work/"prior/hybrid"), ("coupled", work/"coupled")):
                u, w = load(folder, n)
                # Use the raw 3-D samples for both discrete operators.
                with np.load(folder/f"fields_{n}.npz") as data:
                    fields.append((name, data["psi"], data["W_cartesian"], quartic_samples(w, 1)))
                    np.testing.assert_allclose(data["W_cartesian"], longitude(w, azimuth), atol=2e-15, rtol=0)
                    np.testing.assert_allclose(data["psi"], np.broadcast_to(u[..., None], audit.shape), atol=2e-15, rtol=0)
            for factor in (1, 4):
                u, samples = reference_fields(archive, work/"prior/levels", n, factor, hashes)
                fields.append((f"reference_{factor}N", np.broadcast_to(u[..., None], audit.shape),
                               longitude(samples[:, ::8, ::factor], azimuth), samples))
            for field, u, w, samples in fields:
                v1 = momentum(audit, u, w)
                frozen = audit.momentum(u, w, 2)
                for region, mask in masks(audit).items():
                    if region.endswith("_wall"):
                        continue
                    ix = np.flatnonzero(mask[:, 0, 0])
                    independent = continuum(audit, u, samples, n, ix)

                    def rms(value):
                        # Use the same volume weight as the original audit.
                        weight = (audit.r*audit.q*u**6)[ix]
                        return float(np.sqrt(np.sum(weight*np.sum(value**2, axis=0))/weight.sum()))

                    errors = {"v1": v1[:, ix], "frozen": frozen[:, ix], "continuum": independent,
                              "v1_continuum_gap": v1[:, ix]-independent,
                              "frozen_continuum_gap": frozen[:, ix]-independent}
                    for method, residual in errors.items():
                        results.append(dict(n=n, field=field, region=region, method=method,
                                            rms=rms(residual), max_norm=float(np.max(
                                                np.sqrt(np.sum(residual**2, axis=0))))))
    return results


def score(work, rows):
    index = {(r["field"], r["region"], r["n"], r["method"]): r["rms"] for r in rows}
    threshold = 1.5
    recorded = {}
    for name, folder in (("hybrid", "prior/hybrid"), ("coupled", "coupled")):
        report = json.loads((work/folder/"independent-audit.json").read_text())
        assert report["passed"] == 9 and report["total"] == 12
        for row in report["rows"]:
            for region in ("inner", "outer"):
                expected = row["audit"][region+"_boundary_band"]["independent_tensor_M_order2"]["volume_weighted_rms"]
                np.testing.assert_allclose(index[name, region+"_band", row["n"], "frozen"],
                                           expected, atol=2e-14, rtol=0)
        momentum_checks = []
        for check in report["checks"]:
            if check["metric"] != "independent_tensor_M_order2":
                continue
            region = check["band"].replace("_boundary", "")
            coarse, fine = check["pair"]
            order = float(np.log2(index[name, region, coarse, "v1"]/index[name, region, fine, "v1"]))
            momentum_checks.append(dict(band=region, pair=check["pair"], order=order,
                                        coarse_rms=index[name, region, coarse, "v1"],
                                        fine_rms=index[name, region, fine, "v1"],
                                        passed=bool(np.isfinite(order) and order >= threshold)))
        hamiltonian = [v for v in report["checks"] if v["metric"] == "independent_flux_H"]
        recorded[name] = {
            "original_frozen_acceptance": {k: v for k, v in report.items() if k != "rows"},
            "original_frozen_hamiltonian_checks": hamiltonian,
            "new_auditor_diagnostic": {
                "passed": sum(v["passed"] for v in hamiltonian+momentum_checks),
                "total": 12,
                "momentum_checks": momentum_checks,
                "interpretation": ("Threshold arithmetic only. The coupled solver shares the v1 stencil, "
                                   "so its roundoff-scale v1 residual orders are not interpretable as convergence."
                                   if name == "coupled" else "Uses six original frozen Hamiltonian checks.")}}
    return {"auditor": VERSION, "threshold": threshold, "fields": recorded}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path, help="Completed frozen-boundary reproduction directory")
    parser.add_argument("out", type=Path, help="Existing output directory")
    args = parser.parse_args()
    audit_hashes = verified_inputs(args.work)
    prior = ROOT/"research/run43_frozen_boundary/results"
    for checked in prior.iterdir():
        assert digest(checked.read_bytes()) == digest((args.work/"results"/checked.name).read_bytes()), checked
    rows = measure(args.work)
    write_table(args.out/"comparison.csv", rows, ("field", "region", "method"))
    (args.out/"scores.json").write_text(json.dumps(score(args.work, rows), indent=2)+"\n")
    (args.out/"provenance.json").write_text(json.dumps(dict(
        starting_commit=BASE, auditor=VERSION, python=sys.version.split()[0],
        dependencies={name: version(name) for name in ("numpy", "scipy", "pymetis", "mpmath")},
        verified_sha256=audit_hashes,
        verified_prior_results={p.name: digest(p.read_bytes()) for p in sorted(prior.iterdir())}),
        indent=2)+"\n")


if __name__ == "__main__":
    main()
