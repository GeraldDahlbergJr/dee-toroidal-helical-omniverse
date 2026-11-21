"""Strong collocation experiment for the Run 43 coupled constraints."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator, onenormest, splu, spsolve

ROOT = Path(__file__).resolve().parents[2]
for directory in ('research', 'research/run43_frozen_momentum', 'research/run43_hybrid'):
    sys.path.insert(0, str(ROOT / directory))
from physical import periodic
from run43_consistent_momentum import G, Grid, scalar_flux_laplacian
from diagnose import basis, longitude, radial_matrix, RMIN, RMAX, J

GRIDS = (8, 16, 32, 64)


def radial_second(n, h):
    """Centered in every interior row, second-order one-sided at both walls."""
    matrix = sp.diags([1., -2., 1.], [-1, 0, 1], shape=(n+1, n+1), format='lil')
    matrix[0, :4] = [2., -5., 4., -1.]
    matrix[-1, -4:] = [-1., 4., -5., 2.]
    return matrix.tocsr()/h**2


class Momentum:
    def __init__(self, n):
        self.n = n
        h = (RMAX-RMIN)/n
        dt = 2*np.pi/n
        self.r = np.linspace(RMIN, RMAX, n+1)[:, None]
        self.theta = (dt*np.arange(n))[None, :]
        a, da = basis(self.r, self.theta)
        self.size = size = (n+1)*n
        eye = sp.eye(size, format='csr')
        dr = sp.kron(radial_matrix(n+1, h), sp.eye(n), format='csr')
        drr = sp.kron(radial_second(n, h), sp.eye(n), format='csr')
        dtheta = sp.kron(sp.eye(n+1), (periodic(n, 1)-periodic(n, -1))/(2*dt),
                         format='csr')
        spin = np.sin(dt)/dt*J
        first = [sp.kron(np.eye(3), d, format='csr') for d in (dr, dtheta)]
        first.append(sp.kron(spin, eye, format='csr'))

        def diagonal(values):
            return sp.diags(values.ravel(), format='csr')

        def contraction(inner, outer):
            dot = np.sum(inner*outer, axis=0)
            return sp.bmat([
                [diagonal(dot*(i == j) + inner[i]*outer[j] - (2/3)*outer[i]*inner[j])
                 for j in range(3)] for i in range(3)
            ], format='csr')

        self.operator = sp.csr_matrix((3*size, 3*size))
        for i in range(3):
            for j in range(3):
                second = sp.kron(np.eye(3), drr, format='csr') if i == j == 0 else first[i]@first[j]
                self.operator += contraction(a[j], a[i])@second
                self.operator += contraction(da[i][j], a[i])@first[j]

        # Strain used in the Hamiltonian uses the same first derivatives.
        gradients = [
            sum(sp.kron(np.eye(3), diagonal(a[j][i]), format='csr')@first[j]
                for j in range(3)) for i in range(3)
        ]
        div = sum(gradients[i][i*size:(i+1)*size] for i in range(3))
        self.strain = sp.vstack([
            gradients[j][i*size:(i+1)*size] + gradients[i][j*size:(j+1)*size]
            - (2/3)*(i == j)*div for i in range(3) for j in range(3)
        ], format='csr')
        self.free = np.flatnonzero(np.tile((np.arange(size) >= n) & (np.arange(size) < size-n), 3))
        self.fixed = np.setdiff1d(np.arange(3*size), self.free)
        self.matrix = self.operator[self.free][:, self.free].tocsc()
        self.lu = splu(self.matrix)

    def apply(self, w):
        return (self.operator@w.ravel()).reshape(w.shape)

    def solve(self, force, walls=None):
        values = np.zeros(3*self.size)
        rhs = force.ravel()[self.free].copy()
        if walls is not None:
            values[self.fixed] = walls.ravel()[self.fixed]
            rhs -= self.operator[self.free][:, self.fixed]@values[self.fixed]
        x = self.lu.solve(rhs)
        values[self.free] = x
        residual = self.matrix@x-rhs
        return values.reshape(3, self.n+1, self.n), dict(
            momentum_linf=float(np.max(abs(residual))),
            momentum_l2=float(np.linalg.norm(residual)),
            rhs_l2=float(np.linalg.norm(rhs)),
        )

    def stability(self):
        shape = self.matrix.shape
        inverse = LinearOperator(shape, matvec=self.lu.solve,
                                 rmatvec=lambda x: self.lu.solve(x, trans='T'))
        norm_inv = onenormest(inverse, t=1)
        norm_matrix = np.asarray(abs(self.matrix).sum(axis=0)).max()
        return dict(inverse_one_norm_estimate=float(norm_inv),
                    condition_one_norm_estimate=float(norm_matrix*norm_inv))


def coupled(momentum):
    n = momentum.n
    g = Grid(n+1, n, 4, build_vector=False)
    nodes = slice(None, None, 4)
    size = momentum.size
    grid = SimpleNamespace(shape=(n+1, n, 1), size=size, r=g.r[nodes],
                           q=g.q[nodes], t=g.t[nodes])
    lap = scalar_flux_laplacian(grid)
    free = np.arange(n, size-n)
    lap_free = lap[free][:, free]
    psi = np.ones(size)
    spatial, source = g.spatial[nodes], g.S[:, nodes]
    history = []

    def vector():
        w, info = momentum.solve(8*np.pi*G*source*psi**6)
        a2 = np.sum((momentum.strain@w.ravel()).reshape(9, size)**2, axis=0)
        return w, a2, info

    def scalar(a2):
        energy = .5*g.temporal + .5*spatial*psi**-4 + g.U
        residual = lap@(psi-1) + a2*psi**-7/8 + 2*np.pi*G*energy*psi**5
        return residual[free], energy

    w, a2, info = vector()
    for _ in range(30):
        old_psi = psi.copy()
        for steps in range(12):
            residual, energy = scalar(a2)
            if np.max(abs(residual)) < 1e-10:
                break
            derivative = -7*a2*psi**-8/8 + 2*np.pi*G*(-2*spatial+5*energy*psi**4)
            step = spsolve(lap_free+sp.diags(derivative[free]), -residual)
            alpha = 1.
            while np.min(psi[free]+alpha*step) <= 0:
                alpha *= .5
            psi[free] += alpha*step
        else:
            raise RuntimeError('scalar Newton failed')
        updated, a2, info = vector()
        hmax = float(np.max(abs(scalar(a2)[0])))
        change = float(np.max(abs(w-updated)))
        history.append(dict(info, H_linf=hmax, W_change=change, newton_steps=steps,
                            psi_change=float(np.max(abs(psi-old_psi)))))
        w = updated
        if change < 1e-12 and hmax < 1e-9:
            return psi.reshape(n+1, n), w, history
    raise RuntimeError('coupled iteration failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=False)
    for n in GRIDS:
        momentum = Momentum(n)
        psi, w, history = coupled(momentum)
        field = longitude(w, 2*np.pi*np.arange(n)/n)
        path = args.directory/f'fields_{n}.npz'
        np.savez_compressed(path, psi=np.broadcast_to(psi[..., None], field.shape[1:]),
                            W_cartesian=field)
        record = dict(n=n, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      history=history, stability=momentum.stability())
        (args.directory/f'result_{n}.json').write_text(json.dumps(record, indent=2)+'\n')
        print(n, history[-1], flush=True)


if __name__ == '__main__':
    main()
