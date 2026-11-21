"""Versioned strong-form Cartesian momentum audit for uniform toroidal grids."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "reproduce"))
from diagnose_toroidal_boundary_convergence import Audit, G, R, RMIN, RMAX

VERSION = "run43-continuum-momentum-v1"


def radial_first(w, h):
    result = np.empty_like(w)
    result[:, 1:-1] = (w[:, 2:]-w[:, :-2])/(2*h)
    result[:, 0] = (4*(w[:, 1]-w[:, 0])-(w[:, 2]-w[:, 0]))/(2*h)
    result[:, -1] = (4*(w[:, -1]-w[:, -2])-(w[:, -1]-w[:, -3]))/(2*h)
    return result


def radial_second(w, h):
    result = np.empty_like(w)
    result[:, 1:-1] = ((w[:, 2:]-w[:, 1:-1])+(w[:, :-2]-w[:, 1:-1]))/h**2
    for wall, step in ((0, 1), (-1, -1)):
        result[:, wall] = (-5*(w[:, wall+step]-w[:, wall])
                           +4*(w[:, wall+2*step]-w[:, wall])
                           -(w[:, wall+3*step]-w[:, wall]))/h**2
    return result


def metric_vectors(audit):
    er, et, ep = audit.er, audit.et, audit.ep
    r, q = audit.r, audit.q
    c, s = np.cos(audit.t), np.sin(audit.t)
    vectors = er, et/r, ep/q
    # First index is the differentiation direction. r and theta do not
    # rotate e_phi, and q is independent of phi.
    derivatives = ((None, -et/r**2, -c*ep/q**2),
                   (et, -er/r, r*s*ep/q**2),
                   (c*ep, -s*ep/r, (s*et-c*er)/q))
    return vectors, derivatives


def contraction(v, a, b):
    """Apply the symmetric trace-free outer product L(v,a) to b."""
    return (v*np.sum(a*b, axis=0) + a*np.sum(v*b, axis=0)
            -(2/3)*b*np.sum(v*a, axis=0))


def divergence(audit, w):
    """Discretize div(LW) on every radial row for general 3-D Cartesian W."""
    w = np.asarray(w, dtype=np.float64)
    if w.shape != (3, *audit.shape) or min(audit.shape) < 4:
        raise ValueError("expected Cartesian W on a grid with at least four nodes per axis")
    a, da = metric_vectors(audit)

    def derivative(i, value):
        if i == 0:
            return radial_first(value, audit.dr)
        spacing = (audit.dt, audit.dp)[i-1]
        return audit.angular_derivative(value, i+1, spacing, 2)

    first = [derivative(i, w) for i in range(3)]
    result = np.zeros_like(w)
    for i in range(3):
        for j in range(i, 3):
            second = radial_second(w, audit.dr) if i == j == 0 else derivative(j, first[i])
            result += contraction(second, a[j], a[i])
            if i != j:
                # Tensor-product differences commute. Keep both contractions.
                result += contraction(second, a[i], a[j])
                result += contraction(first[i], da[j][i], a[j])
            if da[i][j] is not None:
                result += contraction(first[j], da[i][j], a[i])
    return result


def momentum(audit, psi, w):
    """Return the momentum residual with original source and conformal scaling."""
    psi = np.asarray(psi, dtype=np.float64)
    if psi.shape != audit.shape or not np.isfinite(psi).all() or np.any(psi <= 0):
        raise ValueError("expected finite positive psi on the audit grid")
    if not np.isfinite(w).all():
        raise ValueError("expected finite Cartesian W")
    return psi**-10*divergence(audit, w)-8*np.pi*G*psi**-4*audit.S


def masks(audit):
    r = audit.r[:, 0, 0]
    inner = (r <= .10+1e-12) & (r > RMIN) & (r < RMAX)
    outer = (r >= .20-1e-12) & (r > RMIN) & (r < RMAX)
    index = np.arange(len(r))
    regions = {"inner_wall": index == 0, "inner_first": index == 1,
               "inner_remaining": inner & (index != 1), "inner_band": inner,
               "outer_wall": index == len(r)-1, "outer_first": index == len(r)-2,
               "outer_remaining": outer & (index != len(r)-2), "outer_band": outer}
    return {name: np.broadcast_to(keep[:, None, None], audit.shape)
            for name, keep in regions.items() if keep.any()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field", type=Path, help="NPZ with psi and W_cartesian")
    args = parser.parse_args()
    with np.load(args.field) as data:
        psi, w = data["psi"], data["W_cartesian"]
    audit = Audit(psi.shape)
    residual = momentum(audit, psi, w)
    print(json.dumps({"auditor": VERSION, "grid_shape": list(psi.shape),
                      "momentum": {name: audit.norms(residual, psi, mask)
                                   for name, mask in masks(audit).items()}}, indent=2))


if __name__ == "__main__":
    main()
