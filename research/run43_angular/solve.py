"""Solve the archived Galerkin problem on the two refined angular grids."""
import argparse
import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pymetis
from scipy.sparse.linalg import splu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run43_reference"))
import reference
from compare import digest

GRIDS = (8, 16, 32, 64)


class OrderedLU:
    """Apply a symmetric graph ordering before SuperLU to limit fill."""
    def __init__(self, matrix):
        graph = matrix.tocsr(copy=True)
        graph.setdiag(0)
        graph.eliminate_zeros()
        graph = graph.astype(bool).maximum(graph.T.astype(bool))
        graph.sort_indices()
        permutation, _ = pymetis.nested_dissection(
            pymetis.CSRAdjacency(graph.indptr, graph.indices))
        self.permutation = np.asarray(permutation)
        del graph
        self.factor = splu(matrix[self.permutation, :][:, self.permutation].tocsc(),
                           permc_spec="NATURAL", diag_pivot_thresh=0.0,
                           options={"SymmetricMode": True})

    def solve(self, rhs):
        result = np.empty_like(rhs)
        result[self.permutation] = self.factor.solve(rhs[self.permutation])
        return result


def solve_grid(directory, angles, radial):
    print(f"angles={angles}, radial={radial}", flush=True)
    with patch.object(reference, "splu", OrderedLU), patch.object(
        reference, "spsolve", lambda matrix, rhs: OrderedLU(matrix).solve(rhs)
    ):
        solver = reference.Reference(angles, radial)
        mms = solver.manufactured()
        psi, displacement, history = solver.coupled()
    assert np.isfinite(psi).all() and np.isfinite(displacement).all()
    assert history[-1]["vector_l2"] < 1e-13
    assert history[-1]["scalar_linf"] < 1e-14
    dest = directory / f"a{angles}-r{radial}"
    dest.mkdir(parents=True)
    path = dest / "reduced.npz"
    np.savez_compressed(path, psi=psi, W=displacement)
    record = dict(angles=angles, radial=radial, history=history,
                  manufactured_max_node_error=mms["max_node_error"],
                  manufactured_linear_residual=mms["linear_residual"],
                  sha256=digest(path))
    (dest / "result.json").write_text(json.dumps(record, indent=2) + "\n")
    print(record, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--grids", type=int, nargs="+", choices=GRIDS, default=GRIDS)
    args = parser.parse_args()
    grids = dict.fromkeys((angular*n, radial*n) for n in args.grids
                          for angular in (2, 4) for radial in (4, 8))
    for angles, radial in grids:
        solve_grid(args.directory, angles, radial)
