"""Apply the frozen momentum operator on the original Run 43 audit nodes."""
import argparse
import hashlib
import json
import platform
import subprocess
import sys
from importlib.metadata import version
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import scipy

from solve import GRIDS
from compare import ROOT, Audit, csv_table, digest

ARCHIVE = Path("research/run43_reference/results/Run43_Radial_Reference_Fields.zip")
ARCHIVE_SHA256 = "8dcd1904f098b1e0748078f8fea63ce009c253dcb7f06b75be70004f3f97b4fa"
AUDITOR = Path("reproduce/diagnose_toroidal_boundary_convergence.py")
AUDITOR_BLOB = "5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22"
TERMS = ("reference_sq", "angular_sq", "hybrid_minus_archived_sq",
         "cross_ref_angular", "cross_ref_hybrid", "cross_angular_hybrid")


def restrict(path, n, angles, radial):
    assert angles % n == radial % n == 0
    with np.load(path) as data:
        scalar, vector = data["psi"], data["W"]
    assert scalar.shape == (radial+1, angles)
    assert vector.shape == (3, radial+1, angles)
    u = scalar[::radial//n, ::angles//n]
    wr, wp, wz = vector[:, ::radial//n, ::angles//n, None]
    # These phi nodes coincide with every (angles // n)th solve node.
    phi = np.arange(n)*2*np.pi/n
    c, s = np.cos(phi), np.sin(phi)
    w = np.array(np.broadcast_arrays(wr*c-wp*s, wr*s+wp*c, wz))
    psi = np.broadcast_to(u[:, :, None], (n+1, n, n))
    assert np.isfinite(psi).all() and np.isfinite(w).all() and np.min(psi) > 0
    assert np.all(psi[[0, -1]] == 1) and np.count_nonzero(w[:, [0, -1]]) == 0
    return psi, w


def summarize(weight, h, archived, fine, coarse, fine_half, coarse_half, archived_half):
    def inner(x, y):
        return float(np.sum(weight*np.sum(x*y, axis=0))/weight.sum())

    def norm(x):
        return np.sqrt(inner(x, x))

    def terms(ref, old):
        a, b, c = ref, old-ref, h-old
        np.testing.assert_allclose(a+b+c, h, atol=2e-18, rtol=2e-15)
        return dict(zip(TERMS, (inner(a, a), inner(b, b), inner(c, c),
                                2*inner(a, b), 2*inner(a, c), 2*inner(b, c))))

    result = terms(fine, archived)
    hybrid_sq = inner(h, h)
    closure = hybrid_sq-sum(result.values())
    assert abs(closure) <= 2e-14*hybrid_sq
    result.update(hybrid_sq=hybrid_sq, closure=closure,
                  archived_rms=norm(archived), refined_2n_rms=norm(coarse),
                  refined_4n_rms=norm(fine), angular_2n_to_4n_rms=norm(fine-coarse),
                  angular_2n_to_4n_rms_change=norm(fine)-norm(coarse),
                  archived_to_4n_rms=norm(archived-fine),
                  radial_n_rms=norm(archived-archived_half),
                  radial_2n_rms=norm(coarse-coarse_half),
                  radial_4n_rms=norm(fine-fine_half))
    variants = {"2n": terms(coarse, archived),
                "radial": terms(fine_half, archived_half),
                "radial_4n": terms(fine_half, archived),
                "radial_n": terms(fine, archived_half)}
    return result, variants


def analyze(data, hybrid, out):
    out.mkdir(parents=True)
    frozen = ROOT / "research/run43_reference/results"
    provenance = json.loads((frozen / "provenance.json").read_text())
    hashes = provenance["sha256"]
    for path, expected in hashes.items():
        assert digest(ROOT / path) == expected, path
    assert digest(ROOT / ARCHIVE) == ARCHIVE_SHA256
    assert subprocess.check_output(["git", "hash-object", ROOT / AUDITOR],
                                   text=True).strip() == AUDITOR_BLOB
    previous = json.loads((frozen / "histories.json").read_text())
    archived_hashes = {(row["n"], row["nr"]): row["sha256"] for row in previous}
    existing = json.loads((ROOT / "research/run43_hybrid/measurements.json").read_text())
    rows, histories, pairs = [], {}, []
    by_region = {"first_row": [], "inner_band": []}
    with ZipFile(ROOT / ARCHIVE) as archive:
        for n, expected, old in zip(GRIDS, provenance["hybrid_sha256"],
                                    existing["experiments"]["hybrid"]["levels"]):
            path = hybrid / f"fields_{n}.npz"
            assert digest(path) == expected == old["sha256"]
            with np.load(path) as data_h:
                psi, w = data_h["psi"], data_h["W_cartesian"]
            audit = Audit(psi.shape)
            radius = audit.r[:, 0, 0]
            indices = np.flatnonzero((radius <= .10+1e-12) & (np.arange(n+1) > 0))
            weight = (audit.r*audit.q*psi**6)[indices]
            h = audit.momentum(psi, w, 2)[:, indices]
            residuals = {}
            for angles in (n, 2*n, 4*n):
                for radial in (4*n, 8*n):
                    if angles == n:
                        raw = archive.read(f"n{n}-r{radial}/reduced.npz")
                        assert hashlib.sha256(raw).hexdigest() == archived_hashes[n, radial]
                        fields = restrict(BytesIO(raw), n, angles, radial)
                    else:
                        folder = data / f"a{angles}-r{radial}"
                        info = json.loads((folder / "result.json").read_text())
                        assert (info["angles"], info["radial"]) == (angles, radial)
                        assert digest(folder / "reduced.npz") == info["sha256"]
                        histories[angles, radial] = info
                        fields = restrict(folder / "reduced.npz", n, angles, radial)
                    residuals[angles, radial] = audit.momentum(*fields, 2)[:, indices]
            for region, selection in (("first_row", slice(0, 1)), ("inner_band", slice(None))):
                fixed, half = 8*n, 4*n
                args = [residuals[angles, radial][:, selection] for angles, radial in
                        ((n, fixed), (4*n, fixed), (2*n, fixed),
                         (4*n, half), (2*n, half), (n, half))]
                result, variants = summarize(weight[selection], h[:, selection], *args)
                row = dict(n=n, region=region, **result)
                rows.append(row)
                by_region[region].append((row, variants))
            assert abs(np.sqrt(rows[-1]["hybrid_sq"]) -
                       old["band_rms"]["inner_boundary_band"]["independent_tensor_M_order2"]) < 1e-15

    for region, grid_rows in by_region.items():
        for (a, va), (b, vb) in zip(grid_rows, grid_rows[1:]):
            hybrid_deficit = b["hybrid_sq"]-a["hybrid_sq"]/8
            for term in TERMS:
                deficit = b[term]-a[term]/8
                pairs.append(dict(region=region, coarse=a["n"], fine=b["n"], term=term,
                                  hybrid_order=np.log2(a["hybrid_sq"]/b["hybrid_sq"])/2,
                                  hybrid_gate_deficit=hybrid_deficit, term_gate_deficit=deficit,
                                  signed_share=deficit/hybrid_deficit,
                                  **{f"{key}_deficit_shift": (vb[key][term]-b[term]) -
                                     (va[key][term]-a[term])/8 for key in va}))
            assert abs(sum(p["term_gate_deficit"] for p in pairs[-len(TERMS):]) -
                       hybrid_deficit) < 2e-14*(b["hybrid_sq"]+a["hybrid_sq"]/8)
    csv_table(out / "decomposition.csv", rows)
    csv_table(out / "pairs.csv", pairs)
    (out / "histories.json").write_text(json.dumps(list(histories.values()), indent=2)+"\n")
    completed = subprocess.run([sys.executable, ROOT / "research/audit_run43_arrays.py", hybrid],
                               cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 1 and not completed.stderr, completed.stderr
    report = json.loads((hybrid / "independent-audit.json").read_text())
    assert report["passed"] == 9 and report["total"] == 12 and len(report["checks"]) == 12
    report.pop("rows")
    (out / "hybrid-audit.json").write_text(json.dumps(report, indent=2)+"\n")
    files = [Path(__file__), Path(__file__).with_name("solve.py"), ROOT / ARCHIVE]
    hashes.update({str(path.resolve().relative_to(ROOT)): digest(path) for path in files})
    metadata = dict(base_commit="c29162b5014de0a6ffdc33db5cc0be1304af7bfc",
                    python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                    pymetis=version("pymetis"), auditor_git_blob=AUDITOR_BLOB,
                    sha256=hashes, verified_hybrid_sha256=dict(zip(GRIDS, provenance["hybrid_sha256"])),
                    verified_radial_sha256={f"n{n}-r{r}": archived_hashes[n, r]
                                            for n in GRIDS for r in (4*n, 8*n)})
    (out / "provenance.json").write_text(json.dumps(metadata, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("data", "hybrid", "out"):
        parser.add_argument(name, type=lambda path: Path(path).resolve())
    args = parser.parse_args()
    analyze(args.data, args.hybrid, args.out)
