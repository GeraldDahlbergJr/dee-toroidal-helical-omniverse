"""Extend the archived angular solve with the same fixed radial mesh."""
import argparse
import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run43_angular"))
from solve import GRIDS, OrderedLU, digest, reference, solve_grid


def run(directory, n, factor, radial_factor=8):
    angles, radial = factor*n, radial_factor*n
    if angles <= radial:
        solve_grid(directory, angles, radial)
        return
    # Archived manufactured reporting assumes radial >= angles for its
    # sampling stride. The production iteration has no such restriction.
    print(f"angles={angles}, radial={radial}", flush=True)
    with patch.object(reference, "splu", OrderedLU), patch.object(
        reference, "spsolve", lambda matrix, rhs: OrderedLU(matrix).solve(rhs)
    ):
        solver = reference.Reference(angles, radial)
        psi, displacement, history = solver.coupled()
    assert np.isfinite(psi).all() and np.isfinite(displacement).all()
    assert history[-1]["vector_l2"] < 1e-13
    assert history[-1]["scalar_linf"] < 1e-14
    folder = directory / f"a{angles}-r{radial}"
    folder.mkdir(parents=True)
    path = folder / "reduced.npz"
    np.savez_compressed(path, psi=psi, W=displacement)
    record = dict(angles=angles, radial=radial, history=history, sha256=digest(path))
    (folder / "result.json").write_text(json.dumps(record, indent=2) + "\n")
    print(record, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--grids", type=int, nargs="+", choices=GRIDS, default=GRIDS)
    parser.add_argument("--factors", type=int, nargs="+", default=[2, 4, 8])
    parser.add_argument("--radial-factor", type=int, choices=(4, 8), default=8)
    args = parser.parse_args()
    for n in args.grids:
        for factor in args.factors:
            run(args.directory, n, factor, args.radial_factor)
