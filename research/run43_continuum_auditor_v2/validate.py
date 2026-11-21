"""Run the two targeted Run 43 auditor qualification experiments and archive checks."""
import argparse
import csv
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import platform
import subprocess
from zipfile import ZipFile

import numpy as np

from diagnostics import (EPS, SAFETY_FACTOR, angular_tail, classify, divergence,
                         operator_scale, physical_rhs, rms)
from auditor import Audit, G, ROOT, masks, divergence as v1, momentum

RESEARCH = ROOT / "research"
V1 = RESEARCH / "run43_continuum_auditor_v1/results"
BASE = "59823a57cf1d83c5a792f874c79e891dc6793167"


def pin(path, expected):
    observed = sha256(path.read_bytes()).hexdigest()
    if observed != expected:
        raise ValueError(f"SHA-256 mismatch: {path}")
    return observed


def base_blob(path):
    relative = path.relative_to(ROOT)
    expected = subprocess.check_output(["git", "rev-parse", f"{BASE}:{relative}"], cwd=ROOT)
    observed = subprocess.check_output(["git", "hash-object", path], cwd=ROOT)
    if observed != expected:
        raise ValueError(f"Starting-commit blob mismatch: {path}")
    return observed.decode().strip()


def angular_tests():
    rows, walls, orders = [], [], []
    for n in (16, 32, 64):
        for direction in (1, 2):
            shape = [5, 8, 8]
            shape[direction] = n
            audit = Audit(shape)
            theta = 2*np.pi*np.arange(n)/n
            layout = [1, 1, 1]
            layout[direction] = n
            psi = np.ones(audit.shape)
            mask = np.ones(audit.shape, dtype=bool)
            other = np.zeros((3, *audit.shape))
            other[0] = 1.0
            a = (audit.et/audit.r, audit.ep/audit.q)[direction-1]
            # Hessian of the angle in Cartesian coordinates. A symmetric
            # outer product is contracted directly with the polarization e_x.
            er, et, ep = audit.er, audit.et, audit.ep
            cylindrical_r = np.cos(audit.t)*er - np.sin(audit.t)*et
            if direction == 1:
                coefficient = -np.sin(audit.t)/(audit.r*audit.q)
                lap_angle = coefficient
                hess_x = -(er*et[0] + et*er[0])/audit.r**2 + coefficient*ep*ep[0]
            else:
                lap_angle = 0
                hess_x = -(cylindrical_r*ep[0] + ep*cylindrical_r[0])/audit.q**2
            for category, k in (("low", 1), ("intermediate", n//4),
                                ("near-nyquist", n//2-1), ("nyquist", n//2)):
                wave = (-1.0)**np.arange(n) if k == n//2 else np.cos(k*theta)
                slope = np.zeros(n) if k == n//2 else -k*np.sin(k*theta)
                w = other * wave.reshape(layout)
                exact = (-k*k*(w*np.sum(a*a, axis=0) + a*a[0]*w[0]/3)
                         + slope.reshape(layout)*(other*lap_angle + hess_x/3))
                fd = v1(audit, w)
                fft = divergence(audit, w, spectral=True)
                target = rms(audit, psi, exact, mask)
                row = dict(n=n, direction=("theta", "phi")[direction-1],
                           category=category, mode=k,
                           first_symbol_ratio=float(np.sinc(2*k/n)),
                           fd_to_exact_rms=rms(audit, psi, fd, mask)/target,
                           fd_relative_error=rms(audit, psi, fd-exact, mask)/target,
                           spectral_relative_error=rms(audit, psi, fft-exact, mask)/target,
                           angular_screen=angular_tail(w, direction+1))
                assert row["spectral_relative_error"] < 1e-11, row
                assert (row["angular_screen"]["classification"] == "screen-clear-on-samples") == (
                    category == "low")
                if k == n//2:
                    assert not np.any(fd)
                rows.append(row)
                if category == "nyquist":
                    x = (audit.r-audit.r[0])/(audit.r[-1]-audit.r[0])
                    envelope = other * (x*(1-x))
                    checker = envelope * wave.reshape(layout)
                    interior = mask.copy()
                    interior[[0, -1]] = False
                    centered = v1(audit, checker)
                    spectral = divergence(audit, checker, spectral=True)
                    parity_control = v1(audit, envelope) * wave.reshape(layout)
                    radial_gap = rms(audit, psi, centered-parity_control, interior)
                    assert radial_gap == 0
                    walls.append(dict(n=n, direction=row["direction"],
                                      fd_rms=rms(audit, psi, centered, interior),
                                      spectral_rms=rms(audit, psi, spectral, interior),
                                      fd_gap_to_radial_only_rms=radial_gap))
    for direction in ("theta", "phi"):
        low = [row for row in rows if row["direction"] == direction and row["category"] == "low"]
        for coarse, fine in zip(low, low[1:]):
            order = np.log2(coarse["fd_relative_error"]/fine["fd_relative_error"])
            assert 1.9 < order < 2.1
            orders.append(dict(direction=direction, pair=[coarse["n"], fine["n"]], order=order))
    # A sine at coarse Nyquist has exactly zero nodal values. A half-cell
    # shift gives alternating +/-1, which is measurable even without a solve.
    return {"cosine_modes": rows, "low_mode_error_orders": orders,
            "homogeneous_wall_nyquist": walls, "nyquist_sine_coarse_sample_max": 0.0,
            "nyquist_sine_half_cell_sample_max": float(
                np.max(abs(np.sin(np.pi*(np.arange(64)+0.5)))))}


def roundoff_metrics(audit, psi, w, residual, wide, rhs, regions):
    scale = operator_scale(audit, psi, w, rhs)
    rows = {}
    for region, mask in regions.items():
        value = rms(audit, psi, residual, mask)
        scale_rms = rms(audit, psi, scale, mask)
        floor = SAFETY_FACTOR * EPS * scale_rms
        evaluation_error = rms(audit, psi, np.asarray(residual-wide, dtype=float), mask)
        rows[region] = dict(residual_rms=value, extended_residual_rms=rms(
                                audit, psi, np.asarray(wide, dtype=float), mask),
                            operator_scale_rms=scale_rms,
                            float64_vs_extended_rms=evaluation_error,
                            evaluation_error_eps_scale_ratio=evaluation_error / (EPS * scale_rms),
                            floor_rms=floor, classification=classify(value, floor))
    return rows


def calibrated_tests():
    rows = []
    previous = {}
    for n in (16, 32, 64):
        audit = Audit((n+1, n, 8))
        for amplitude in (2.0**-20, 1.0, 2.0**20):
            w = amplitude * audit.er * (1 + audit.r) * (1 + np.sin(3*audit.t)/8)
            psi = np.broadcast_to(1 + audit.r/8, audit.shape)
            lhs = psi**-10 * v1(audit, w)
            lhs_wide = psi.astype(np.longdouble)**-10 * divergence(audit, w.astype(np.longdouble))
            mask = masks(audit)["inner_band"]
            scale = rms(audit, psi, operator_scale(audit, psi, w, lhs), mask)
            # Prescribed RHS includes a halving-squared signal or zero signal.
            # The nonzero signal has analytic order 2 after removing input scaling.
            signals = {"roundoff-limited": 0.0,
                       "near-floor": 4 * SAFETY_FACTOR * EPS * scale,
                       "above-floor": amplitude * 2.0**-16 * (16/n)**2}
            parity = (-1.0)**np.arange(len(audit.r))[None, :, None, None]
            adjacent = np.nextafter(w, np.where(parity > 0, np.inf, -np.inf))
            input_change = rms(audit, psi, psi**-10 * v1(audit, adjacent-w), mask)
            for expected, signal in signals.items():
                exact_rhs = lhs_wide.copy()
                exact_rhs[0] -= signal
                rhs = exact_rhs.astype(float)
                metrics = roundoff_metrics(audit, psi, w, lhs-rhs, lhs_wide-exact_rhs,
                                           rhs, {"inner_band": mask})["inner_band"]
                assert metrics["classification"] == expected, (n, amplitude, metrics)
                order = None
                if expected == "above-floor":
                    value = metrics["residual_rms"]
                    if amplitude in previous:
                        order = float(np.log2(previous[amplitude]/value))
                        assert abs(order-2) < 1e-4
                    previous[amplitude] = value
                rows.append(dict(n=n, amplitude=amplitude, experiment=expected, signal=signal,
                                 one_ulp_input_change_rms=input_change, order=order, **metrics))
    return rows


def longitude(w, phi):
    x, y, z = w[..., None]
    return np.stack((x*np.cos(phi)-y*np.sin(phi), x*np.sin(phi)+y*np.cos(phi),
                     z + np.zeros_like(phi)))


def inspect_field(audit, psi, w):
    rhs = physical_rhs(audit, psi)
    residual = momentum(audit, psi, w)
    u = psi.astype(np.longdouble)
    # Hold the sampled geometry and physical source S fixed. Compute the
    # conformal factors, stencil and contractions in extended precision.
    wide = u**-10 * divergence(audit, w.astype(np.longdouble)) - (
        8 * np.longdouble(np.pi) * G * u**-4 * audit.S)
    regions = {key: mask for key, mask in masks(audit).items() if not key.endswith("_wall")}
    for key, row in roundoff_metrics(audit, psi, w, residual, wide, rhs, regions).items():
        yield key, row


def archives(hybrid, refined):
    recorded = json.loads((V1/"provenance.json").read_text())["verified_sha256"]
    checked = {str(path.relative_to(ROOT)): base_blob(path) for path in (
        V1/"comparison.csv", V1/"scores.json", V1/"provenance.json",
        V1.parent/"auditor.py", V1.parent/"REPORT.md")}
    frozen = "reproduce/diagnose_toroidal_boundary_convergence.py"
    path = ROOT/frozen
    checked[frozen] = pin(path, recorded[frozen])
    assert base_blob(path) == "5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22"
    comparison, reported_orders = {}, {}
    with (V1/"comparison.csv").open() as stream:
        for row in csv.DictReader(stream):
            key = row["field"], int(row["n"]), row["region"]
            comparison.setdefault(key, {})[row["method"]] = float(row["rms"])
            if row["method"] == "continuum" and row["order"]:
                reported_orders[key] = float(row["order"])
    history = RESEARCH/"run43_asymptotic/results/histories.json"
    checked[str(history.relative_to(ROOT))] = pin(history, recorded[str(history.relative_to(ROOT))])
    levels = json.loads(history.read_text())
    member_manifest = RESEARCH/"run43_reference/results/histories.json"
    checked[str(member_manifest.relative_to(ROOT))] = base_blob(member_manifest)
    native = {(r["n"], r["nr"]): r["sha256"] for r in json.loads(member_manifest.read_text())}
    high = {(r["angles"], r["radial"]): r["sha256"] for r in levels if "angles" in r}
    coupled = RESEARCH/"run43_coupled_recovery/results/fields"
    checked["coupled/independent-audit.json"] = pin(
        coupled/"independent-audit.json", recorded["coupled/independent-audit.json"])
    archive_path = RESEARCH/"run43_reference/results/Run43_Radial_Reference_Fields.zip"
    archive_key = str(archive_path.relative_to(ROOT))
    checked[archive_key] = pin(archive_path, recorded[archive_key])
    measures, screens = {}, []
    with ZipFile(archive_path) as archive:
        for n in (8, 16, 32, 64):
            audit = Audit((n+1, n, n))
            for name, directory, key in (("coupled", coupled, "coupled"),
                                          ("hybrid", hybrid, "prior/hybrid")):
                if directory is None:
                    continue
                filename = f"fields_{n}.npz"
                checked[f"{key}/{filename}"] = pin(directory/filename, recorded[f"{key}/{filename}"])
                with np.load(directory/filename) as data:
                    psi, w = data["psi"], data["W_cartesian"]
                screens.append(screen(name, n, psi, w))
                for region, row in inspect_field(audit, psi, w):
                    measures[name, n, region] = row
            for factor in (1, 4):
                if factor == 1:
                    member = f"n{n}-r{8*n}/reduced.npz"
                    raw = archive.read(member)
                    expected = native[n, 8*n]
                elif refined is not None:
                    member = f"a{4*n}-r{8*n}/reduced.npz"
                    raw = (refined/member).read_bytes()
                    expected = high[4*n, 8*n]
                else:
                    continue
                observed = sha256(raw).hexdigest()
                if observed != expected:
                    raise ValueError(f"Reference member mismatch: {member}")
                checked[member] = observed
                name = f"reference_{factor}N"
                with np.load(BytesIO(raw)) as data:
                    u, samples = data["psi"], data["W"]
                psi = np.broadcast_to(u[::8, ::factor, None], audit.shape)
                w = longitude(samples[:, ::8, ::factor], audit.p[0, 0])
                screens.append(screen(name, n, psi, w))
                screens[-1]["native_theta_w"] = angular_tail(samples, 2)
                screens[-1]["native_theta_psi"] = angular_tail(u, 1)
                for region, row in inspect_field(audit, psi, w):
                    measures[name, n, region] = row
    rows = []
    for (field, n, region), data in sorted(comparison.items()):
        result = measures.get((field, n, region))
        if result:
            np.testing.assert_allclose(result["residual_rms"], data["v1"], atol=2e-14, rtol=0)
            if field == "coupled":
                assert result["classification"] == "roundoff-limited"
                assert result["extended_residual_rms"] <= result["floor_rms"]
            if field == "reference_1N":
                assert result["classification"] == "above-floor"
        prior = comparison.get((field, n//2, region))
        prior_check = measures.get((field, n//2, region), {})
        current_check = result or {}
        order = None
        if prior and all(check.get("classification") == "above-floor"
                         for check in (prior_check, current_check)):
            order = float(np.log2(prior["v1"]/data["v1"]))
        rows.append(dict(field=field, n=n, region=region, historical_v1_rms=data["v1"],
                         roundoff_evidence="array-verified" if result else "summary-only",
                         roundoff=result, v1_order=order,
                         independent_residual_rms=data["continuum"],
                         archived_independent_order=reported_orders.get((field, n, region)),
                         v1_continuum_gap_rms=data["v1_continuum_gap"],
                         gap_exceeds_independent_rms=data["v1_continuum_gap"] > data["continuum"],
                         gap_over_independent=data["v1_continuum_gap"]/data["continuum"]))
    return {"checked_hashes": checked, "angular_screens": screens, "regions": rows,
            "historical_acceptance": {name: [field["original_frozen_acceptance"][key]
                                             for key in ("passed", "total")]
                                      for name, field in json.loads(
                                          (V1/"scores.json").read_text())["fields"].items()}}


def screen(name, n, psi, w):
    return dict(field=name, n=n,
                theta_w=angular_tail(w, 2), phi_w=angular_tail(w, 3),
                theta_psi=angular_tail(psi, 1), phi_psi=angular_tail(psi, 2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--hybrid", type=Path, help="Directory containing the original hybrid NPZs")
    parser.add_argument("--refined", type=Path, help="Directory containing original 4N level folders")
    args = parser.parse_args()
    if np.finfo(np.longdouble).eps >= EPS:
        raise RuntimeError("Extended precision must have more significand bits than float64")
    result = dict(base=BASE, python=platform.python_version(), numpy=np.__version__,
                  extended_epsilon=float(np.finfo(np.longdouble).eps),
                  angular=angular_tests(), roundoff=calibrated_tests(),
                  archives=archives(args.hybrid, args.refined))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
