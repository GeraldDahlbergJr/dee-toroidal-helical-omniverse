"""Run 43 physical-source SBP weak-flux candidate (experimental, unpromoted).

This replaces the collocated vector momentum *solve* operator with a
weighted weak strain-flux operator on the actual toroidal shell. The matter
source, nonlinear scalar equation, Dirichlet traces, frozen audit and
acceptance threshold are not changed.

The first-derivative radial closure is SBP(2,1). The full weak operator is
energy symmetric in the weighted scalar product. The vector flux is the
Cartesian trace-free strain sigma=L W, i.e. an auxiliary mixed flux
eliminated locally by sigma=L W.

IMPORTANT: Because homogeneous radial Dirichlet values are eliminated
strongly in the current solver, the symmetric Dirichlet SAT trace terms
vanish in the reduced interior equations. This is a mixed-flux/SBP
candidate, but is NOT a complete weak-SAT enforcement implementation.
Do not claim recovery without the unchanged 12-check audit.
"""
import os
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator
import run43_consistent_momentum as original
from run43_consistent_momentum import TensorGrid
from validate_toroidal_helical_dee_source import RMIN, RMAX


class MixedFluxSBPGrid(TensorGrid):
    def __init__(self, n, stabilization=0.25):
        # Builds the original geometry, source and scalar flux Laplacian.
        super().__init__(n, stabilization=stabilization)
        nr, nt, np_ = self.shape
        dr = (RMAX - RMIN) / (nr - 1)
        dtheta = 2 * np.pi / nt
        dphi = 2 * np.pi / np_

        # Second-order interior / first-order endpoint diagonal-norm SBP D_r.
        radial = sp.lil_matrix((nr, nr))
        radial[0, 0:2] = [-1, 1]
        radial[-1, -2:] = [-1, 1]
        for k in range(1, nr - 1):
            radial[k, k - 1:k + 2] = [-0.5, 0, 0.5]
        radial = (radial.tocsr() / dr)
        hr = np.ones(nr) * dr
        hr[0] *= 0.5
        hr[-1] *= 0.5
        Hrad = sp.diags(hr)
        B = sp.diags([-1] + [0] * (nr - 2) + [1])
        assert np.max(np.abs((Hrad @ radial + radial.T @ Hrad - B).toarray())) < 1e-12

        it = sp.eye(nt, format="csr")
        ip = sp.eye(np_, format="csr")
        self.d[0] = sp.kron(sp.kron(radial, it, format="csr"), ip, format="csr")
        er, et, ep = self.basis
        diag = sp.diags
        self.grad = [
            (diag(er[i]) @ self.d[0]
             + diag(et[i] / self.r) @ self.d[1]
             + diag(ep[i] / self.q) @ self.d[2]).tocsr()
            for i in range(3)
        ]
        weight = np.broadcast_to(
            hr[:, None, None], self.shape).ravel() * dtheta * dphi * self.r * self.q
        assert np.min(weight) > 0
        self.sbp_weight = weight
        self.radial_sbp_defect = float(np.max(
            np.abs((Hrad @ radial + radial.T @ Hrad - B).toarray())))

        def apply(x):
            w = np.zeros((3, self.size))
            w[:, self.idx] = x.reshape(3, -1)
            sigma = self.longitudinal(w)
            # Weak div sigma = -(grad)^*_weighted sigma, using Euclidean
            # Cartesian strain. Strong homogeneous radial boundary values
            # are already eliminated by self.idx / self.inject.
            y = np.asarray([
                -sum(self.grad[j].T @ (weight * sigma[i, j])
                     for j in range(3)) / weight
                for i in range(3)
            ])
            # Preserve the separately documented O(h_angle^2) angular
            # checkerboard damping. It is an added experimental term and
            # must not be confused with the mixed-flux identity.
            if stabilization:
                h = dtheta
                for i in range(3):
                    v = w[i].reshape(self.shape)
                    d4theta = (np.roll(v, 2, axis=1)
                               - 4 * np.roll(v, 1, axis=1)
                               + 6 * v - 4 * np.roll(v, -1, axis=1)
                               + np.roll(v, -2, axis=1)) / h**2
                    hp = dphi
                    d4phi = (np.roll(v, 2, axis=2)
                             - 4 * np.roll(v, 1, axis=2)
                             + 6 * v - 4 * np.roll(v, -1, axis=2)
                             + np.roll(v, -2, axis=2)) / hp**2
                    y[i] -= stabilization * (
                        d4theta.ravel() / self.r**2
                        + d4phi.ravel() / self.q**2)
            return y[:, self.idx].ravel()

        m = 3 * len(self.idx)
        self.vector = LinearOperator((m, m), matvec=apply, dtype=float)


def main():
    original.TensorGrid = MixedFluxSBPGrid
    original.output_dir = lambda: (
        Path(__file__).resolve().parent / "results" /
        os.environ.get("DEE_EXPERIMENT", "mixed-flux-sbp-physical"))
    original.main()


if __name__ == "__main__":
    main()
