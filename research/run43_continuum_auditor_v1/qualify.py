"""Verify v1 against analytic divergence in Cartesian and cylindrical frames."""
import argparse
from pathlib import Path
import sys

import numpy as np

from auditor import Audit, G, R, ROOT, RMIN, RMAX, masks, momentum, radial_second
sys.path.insert(0, str(ROOT/"research/run43_asymptotic"))
from operator_eval import divergence, longitude
sys.path.insert(0, str(ROOT/"research/run43_frozen_momentum"))
from diagnose import write_table


def homogeneous(grid):
    """Equivariant sine envelope with fully analytic cylindrical jets."""
    r, t, p = grid.r[..., 0], grid.t[..., 0], grid.p[0, 0]
    k = np.pi/(RMAX-RMIN)
    f = np.sin(k*(r-RMIN))
    f[[0, -1]] = 0
    fr = k*np.cos(k*(r-RMIN))
    g = np.concatenate((np.cos(2*t), np.sin(t), 1+.2*np.cos(t)))[:, None]
    gt = np.concatenate((-2*np.sin(2*t), np.cos(t), -.2*np.sin(t)))[:, None]
    gtt = np.concatenate((-4*np.cos(2*t), -np.sin(t), -.2*np.cos(t)))[:, None]
    w = f*g
    exact = divergence(w, fr*g, -k*k*w, f*gt, fr*gt, f*gtt, r, t, R)
    return longitude(w, p), longitude(exact, p)


def cartesian_waves(grid):
    """Lifted non-equivariant waves with Cartesian Hessians known exactly."""
    r, t, p = grid.r, grid.t, grid.p
    xyz = (grid.q*np.cos(p), grid.q*np.sin(p),
           np.broadcast_to(r*np.sin(t), grid.shape))
    wavevectors = np.array([[1., 2., 3.], [-2., 3., 1.], [3., 1., -2.]])
    phases = np.einsum("ij,j...->i...", wavevectors, xyz)
    w = np.stack((np.sin(phases[0]), np.cos(phases[1]), np.sin(phases[2])))
    # For each W_j, ∂_i ∂_j W_j = −k_ji k_jj W_j.
    grad_div = -np.einsum("ji,j,j...->i...", wavevectors, wavevectors.diagonal(), w)
    laplacian = -np.sum(wavevectors**2, axis=1)[:, None, None, None]*w
    assert np.linalg.norm(w[:, 0]) > 0 and np.linalg.norm(w[:, -1]) > 0
    return w, laplacian + grad_div/3


def verify_radial_stencil():
    x = np.linspace(0, 1, 17)
    for degree in range(4):
        w = (x**degree)[None, :, None, None]
        expected = degree*(degree-1)*x**max(0, degree-2)
        np.testing.assert_allclose(radial_second(w, x[1])[:, :, 0, 0],
                                   expected[None], atol=1e-13, rtol=0)


def qualify():
    verify_radial_stencil()
    rows = []
    for n in (8, 16, 32, 64, 128):
        grid = Audit((n+1, n, n))
        psi = np.broadcast_to(1+.03*np.cos(grid.t)+.02*np.sin(grid.p), grid.shape)
        for case in (homogeneous, cartesian_waves):
            w, exact = case(grid)
            expected = psi**-10*exact-8*np.pi*G*psi**-4*grid.S
            actual = momentum(grid, psi, w)
            error = actual-expected
            for region, mask in masks(grid).items():
                norm = grid.norms(error, psi, mask)
                rows.append(dict(n=n, case=case.__name__, region=region,
                                 rms=norm["volume_weighted_rms"], max_error=norm["max_vector_norm"]))
    for row in rows:
        if row["n"] == 128:
            prev = next(p for p in rows if p["n"] == 64
                        and (p["case"], p["region"]) == (row["case"], row["region"]))
            for metric in ("rms", "max_error"):
                assert np.log2(prev[metric]/row[metric]) > 1.8, (row, metric, prev)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    write_table(args.out, qualify(), ("case", "region"))
