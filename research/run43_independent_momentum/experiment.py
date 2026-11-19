"""Quadrature Galerkin momentum pilot on archived original-source Run 43 fields."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, gmres, splu

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "research"), str(ROOT / "reproduce")]
from validate_original_physics import Independent
from dee_stress_energy_3p1 import direct_tensor_check
from diagnose_toroidal_boundary_convergence import Audit
from solve_toroidal_helical_coupled_constraints import Grid, G

ARCHIVE = ROOT / "research/results/coordinated-replay-20261103-1831/baseline"
ROTATION = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]])


def q2_axis(n, start, end, periodic, order):
    x, weights = leggauss(order)
    h = (end - start) / n
    n_cells = n // 2
    cols = (2 * np.arange(n_cells)[:, None] + np.arange(3)) % (n + 1 - periodic)
    rows = np.broadcast_to(np.arange(n_cells * order).reshape(n_cells, order, 1),
                           (n_cells, order, 3))
    cols = np.broadcast_to(cols[:, None, :], rows.shape)
    shape = (n_cells * order, n + 1 - periodic)
    values = np.column_stack((x * (x - 1) / 2, 1 - x*x, x * (x + 1) / 2))
    slopes = np.column_stack((x - .5, -2*x, x + .5)) / h
    matrices = [
        sparse.csr_matrix((np.tile(table, (n_cells, 1)).ravel(),
                           (rows.ravel(), cols.ravel())), shape=shape)
        for table in (values, slopes)
    ]
    points = start + h * (2*np.arange(n_cells)[:, None] + 1 + x)
    return *matrices, points.ravel(), np.tile(weights*h, n_cells)


class Galerkin:
    def __init__(self, n, order, mode=0):
        br, dr, r, wr = q2_axis(n, .05, .25, False, order)
        bt, dt, t, wt = q2_axis(n, 0, 2*np.pi, True, order)
        self.n = n
        self.mode = mode
        self.active = np.arange((n+1)*n).reshape(n+1, n)[1:-1].ravel()
        self.basis = sparse.kron(br, bt, format="csr")
        self.B = self.basis[:, self.active]
        self.Dr = sparse.kron(dr, bt, format="csr")[:, self.active]
        self.Dt = sparse.kron(br, dt, format="csr")[:, self.active]
        r, t = (a.ravel() for a in np.meshgrid(r, t, indexing="ij"))
        c, s = np.cos(t), np.sin(t)
        q = 1 + r*c
        self.weights = np.outer(wr, wt).ravel() * r*q
        zero = self.B * 0
        angular = self.B.multiply((1/q)[:, None])
        dx = self.Dr.multiply(c[:, None]) - self.Dt.multiply((s/r)[:, None])
        dz = self.Dr.multiply(s[:, None]) + self.Dt.multiply((c/r)[:, None])
        rotation = ROTATION + 1j*mode*np.eye(3) if mode else ROTATION
        gradients = []
        for i in range(3):
            gradients.append([
                sparse.hstack([dx if i == j else zero for j in range(3)], format="csr"),
                sparse.hstack([angular * rotation[i, j] for j in range(3)], format="csr"),
                sparse.hstack([dz if i == j else zero for j in range(3)], format="csr"),
            ])
        div = sum(gradients[i][i] for i in range(3))
        stiffness = sparse.block_diag([self.gram(dx) + self.gram(dz)]*3, format="csr")
        stiffness += sum(self.gram(gradients[i][1]) for i in range(3))
        self.vector_lap = stiffness.tocsc()
        self.K = (stiffness + self.gram(div)/3).tocsr()
        self.direction = np.array([-s/r, -3/q, c/r])
        self.source = -.030316 * self.direction

    def gram(self, matrix):
        return matrix.conjugate().T @ matrix.multiply(self.weights[:, None])

    def solve(self, psi):
        if self.mode != 0:
            raise ValueError("this load is defined only for the physical axisymmetric mode")
        factor = splu(self.vector_lap)
        preconditioner = LinearOperator(self.K.shape, matvec=factor.solve, dtype=self.K.dtype)
        values = self.basis @ psi.ravel()
        if values.min() <= 0:
            raise ValueError("nonpositive conformal factor at load quadrature")
        rhs = np.concatenate([self.B.T @ (-8*np.pi*G*self.weights*values**6*v)
                              for v in self.source])
        count = []
        x, info = gmres(self.K, rhs, M=preconditioner, rtol=1e-10, atol=1e-13,
                        restart=60, maxiter=100, callback=count.append, callback_type="pr_norm")
        res = self.K @ x - rhs
        rhs_norm = np.linalg.norm(rhs)
        threshold = max(1e-13, 1e-10*rhs_norm)
        if info or np.linalg.norm(res) > threshold:
            raise RuntimeError(f"GMRES info={info}, true residual={np.linalg.norm(res)}")
        fitted = np.zeros((3, self.n+1, self.n))
        fitted.reshape(3, -1)[:, self.active] = x.reshape(3, -1)
        return fitted, dict(iterations=len(count), info=info, rhs_l2=rhs_norm,
                            residual_l2=np.linalg.norm(res), residual_max=np.max(abs(res)),
                            relative=np.linalg.norm(res)/rhs_norm, threshold=threshold)


def cartesian(w):
    n = w.shape[-1]
    phi = np.arange(n)*2*np.pi/n
    c, s = np.cos(phi), np.sin(phi)
    x, y, z = w[..., None]
    return np.array([c*x-s*y, s*x+c*y, np.broadcast_to(z, x.shape[:-1]+(n,))])


def reduce_fields(psi, w):
    n = psi.shape[-1]
    phi = np.arange(n)*2*np.pi/n
    c, s = np.cos(phi), np.sin(phi)
    x, y, z = w
    reduced = np.array([np.mean(c*x+s*y, axis=-1), np.mean(-s*x+c*y, axis=-1),
                        np.mean(z, axis=-1)])
    mean = np.mean(psi, axis=-1)
    error = dict(psi=float(np.max(abs(psi-mean[..., None]))),
                 W=float(np.max(abs(w-cartesian(reduced)))))
    return mean, error


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strong(ind, w):
    """Expanded chain-rule Hessian using an independent radial/Fourier stencil."""
    r, t, q = ind.r, ind.t, ind.q
    c, s = np.cos(t), np.sin(t)
    er, et, ep = ind.er, ind.et, ind.ep
    reciprocal = (er, et/r, ep/q)
    derivatives = (
        (np.zeros_like(er), -et/r**2, -c*ep/q**2),
        (et, -er/r, r*s*ep/q**2),
        (c*ep, -s*ep/r, -(c*er-s*et)/q),
    )
    d = (ind.dr, lambda u: ind.angle(u, 1), lambda u: ind.angle(u, 2))
    d2 = (lambda u: np.einsum("ij,jkl->ikl", ind.d2, u),
          lambda u: ind.angle(u, 1, 2), lambda u: ind.angle(u, 2, 2))
    first = [[differentiate(u) for u in w] for differentiate in d]
    blocks = {name: np.zeros_like(w) for name in ("radial", "angular", "mixed")}
    for a, outer in enumerate(reciprocal):
        for b, inner in enumerate(reciprocal):
            name = "mixed" if (a == 0) != (b == 0) else ("radial" if a == 0 else "angular")
            for j in range(3):
                second = d2[a](w[j]) if a == b else d[a](first[b][j])
                coefficient_gradient = derivatives[a][b]
                blocks[name][j] += (np.sum(outer*inner, axis=0)*second
                                    + np.sum(outer*coefficient_gradient, axis=0)*first[b][j])
                blocks[name] += outer*(inner[j]*second+coefficient_gradient[j]*first[b][j])/3
    return blocks


def archived_operators(n, psi, fields, independent_blocks, masks):
    # Apply the original exact sparse stencils without building their LU.
    grid = Grid(n+1, n, n, build_vector=False)
    rhs = 8*np.pi*G*psi**6*grid.S.reshape((3,) + psi.shape)
    output = {}
    for label, w in fields.items():
        vector = w.reshape(3, -1)
        dw = [[grid.grad[a] @ v for v in vector] for a in range(3)]
        div = sum(dw[i][i] for i in range(3))
        expanded = np.array([grid.lap @ v + grid.grad[i] @ div/3
                             for i, v in enumerate(vector)]).reshape(w.shape)
        tensor = grid.longitudinal(vector)
        flux = np.array([sum(grid.grad[j] @ tensor[i, j] for j in range(3))
                         for i in range(3)]).reshape(w.shape)
        h = 2*np.pi/n
        for axis, scale in ((2, grid.r**2), (3, grid.q**2)):
            fourth = (np.roll(w, 2, axis)-4*np.roll(w, 1, axis)+6*w
                      -4*np.roll(w, -1, axis)+np.roll(w, -2, axis))
            flux -= fourth/(4*h*h*scale.reshape(psi.shape))
        reference = sum(independent_blocks[label].values())
        audit = Audit(psi.shape)
        output[label] = {}
        for name, action in (("original_expanded", expanded), ("archived_stabilized", flux)):
            r = (action-rhs)[:, 1:-1]
            output[label][name] = {
                "residual_l2": np.linalg.norm(r),
                "relative": np.linalg.norm(r)/np.linalg.norm(rhs[:, 1:-1]),
                "difference_from_independent_action": {
                    region: audit.norms(psi**-10*(action-reference), psi, mask)
                    for region, mask in masks.items() if not region.endswith("endpoint")
                },
            }
    return output


def projection_check(ind, psi):
    """Check all original projections against the direct spacetime tensor."""
    direction = ind.et/ind.r - 3*ind.ep/ind.q
    maximum = np.zeros(4)
    n = ind.n
    for i in (0, 1, n//2, n-1, n):
        for j in (0, n//4, n//2):
            k = (i+j) % n
            v = direction[:, i, j, k]
            scale = psi[i, j, k]**4
            gamma = scale*np.eye(3)
            d = np.array([np.zeros(3), .25*v, .20*v])
            direct = direct_tensor_check(np.sqrt(.945), np.array([0, .12, -.04]), d, gamma)
            potential = (.945-1)**2
            half_spatial = .068614375*(v@v)/scale
            pressure = .0127638-2*potential-half_spatial
            stress = 2*.068614375*np.outer(v, v)+pressure*gamma
            analytic = (.0127638+half_spatial, -.030316*v, stress,
                        2*half_spatial+3*pressure)
            maximum = np.maximum(maximum, [np.max(abs(a-b)) for a, b in zip(direct, analytic)])
    return dict(zip(("E_max_abs", "S_i_max_abs", "S_ij_max_abs", "trace_max_abs"), maximum))


def run(out):
    from scipy.sparse.linalg import eigsh

    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    rows = []
    nested = []
    previous = None
    sources = [Path(__file__), ROOT/"research/validate_original_physics.py",
               ROOT/"research/audit_run43_arrays.py", ROOT/"research/run43_consistent_momentum.py",
               ROOT/"reproduce/dee_stress_energy_3p1.py",
               ROOT/"reproduce/validate_toroidal_helical_dee_source.py",
               ROOT/"reproduce/diagnose_toroidal_boundary_convergence.py",
               ROOT/"reproduce/solve_toroidal_helical_coupled_constraints.py",
               ARCHIVE/"executed_solver.py", ARCHIVE/"settings.json"]
    provenance = dict(python=sys.version, platform=platform.platform(), numpy=np.__version__,
                      scipy=scipy.__version__, base=subprocess.check_output(
                          ["git", "rev-parse", "origin/run43-final-certification-20261107"],
                          text=True).strip(),
                      sources={str(p.relative_to(ROOT)): sha(p) for p in sources})
    (out/"environment.json").write_text(json.dumps(provenance, indent=2)+"\n")
    for n in (8, 16, 32, 64):
        path = ARCHIVE/f"fields_{n}.npz"
        digest = sha(path)
        assert digest == json.loads((ARCHIVE/f"result_{n}.json").read_text())["sha256"]
        with np.load(path) as data:
            psi, archived = data["psi"], data["W_cartesian"]
        assert np.array_equal(psi[[0, -1]], np.ones_like(psi[[0, -1]]))
        assert not np.count_nonzero(archived[:, [0, -1]])
        mean, symmetry = reduce_fields(psi, archived)
        galerkin = Galerkin(n, 4)
        w, algebra = galerkin.solve(mean)
        w6, algebra6 = Galerkin(n, 6).solve(mean)
        candidate = cartesian(w)
        fields = {"archive": archived, "candidate": candidate}
        np.savez_compressed(out/f"fields_{n}.npz", psi=psi, W_cartesian=candidate)
        (out/f"result_{n}.json").write_text(json.dumps({"sha256": sha(out/f"fields_{n}.npz")})+"\n")
        v0 = np.random.default_rng(43).normal(size=galerkin.K.shape[0])
        lo = eigsh(galerkin.K, k=1, sigma=0, which="LM", v0=v0,
                   return_eigenvectors=False)[0]
        hi = eigsh(galerkin.K, k=1, which="LA", v0=v0, return_eigenvectors=False)[0]
        ind = Independent(n)
        audit = Audit(psi.shape)
        r = ind.r.ravel()
        planes = np.arange(n+1)
        bands = {
            "interior": (planes > 0) & (planes < n),
            "inner": (planes > 0) & (r <= .10+1e-12),
            "bulk": (r >= .10-1e-12) & (r <= .20+1e-12),
            "outer": (planes < n) & (r >= .20-1e-12),
            "first": planes == 1, "last": planes == n-1,
            "inner_endpoint": planes == 0, "outer_endpoint": planes == n,
        }
        masks = {key: np.broadcast_to(band[:, None, None], psi.shape) for key, band in bands.items()}
        blocks = {key: strong(ind, value) for key, value in fields.items()}
        rhs = 8*np.pi*G*psi**6*(-.030316)*(ind.et/ind.r-3*ind.ep/ind.q)
        row = dict(n=n, archive_sha256=digest, symmetry_defect=symmetry, gmres=algebra,
                   six_point_gmres=algebra6, quadrature_W_max=np.max(abs(w6-w)),
                   direct_spacetime_projections=projection_check(ind, psi),
                   stiffness_eigenvalues=[lo, hi], stiffness_condition=hi/lo,
                   supplementary={key: ind.residuals(psi, value) for key, value in fields.items()},
                   archive_operator_comparison=archived_operators(n, psi, fields, blocks, masks),
                   boundary_and_interior={})
        for name, mask in masks.items():
            row["boundary_and_interior"][name] = {
                **{key: audit.norms(psi**-10*(sum(action.values())-rhs), psi, mask)
                   for key, action in blocks.items()},
                "candidate_minus_archive_W": audit.norms(candidate-archived, psi, mask),
                "action_differences": {
                    key: audit.norms(psi**-10*(blocks["candidate"][key]-blocks["archive"][key]),
                                     psi, mask) for key in blocks["candidate"]
                }
            }
        if previous is not None:
            coarse_psi, coarse_fields, coarse_blocks, coarse_masks = previous
            coarse_audit = Audit(coarse_psi.shape)
            paired = {}
            for label in fields:
                paired[label] = {}
                for region, mask in coarse_masks.items():
                    coarse_quantities = {"W": coarse_fields[label], **coarse_blocks[label]}
                    fine_quantities = {"W": fields[label], **blocks[label]}
                    paired[label][region] = {
                        key: coarse_audit.norms(
                            value-fine_quantities[key][:, ::2, ::2, ::2], coarse_psi, mask)
                        for key, value in coarse_quantities.items()
                    }
            nested.append(dict(pair=[n//2, n], differences=paired))
        previous = psi, fields, blocks, masks
        rows.append(row)
        print(json.dumps(dict(n=n, elapsed_s=time.monotonic()-start, gmres=algebra,
                              condition=hi/lo, quadrature_W_max=row["quadrature_W_max"])),
              flush=True)
        (out/"measurements.json").write_text(json.dumps(rows, indent=2)+"\n")
    (out/"nested.json").write_text(json.dumps(nested, indent=2)+"\n")
    orders = {}
    for region in rows[0]["boundary_and_interior"]:
        orders[region] = {}
        for key in ("archive", "candidate"):
            values = [row["boundary_and_interior"][region][key]["volume_weighted_rms"] for row in rows]
            orders[region][key] = dict(values=values, orders=np.log2(
                np.array(values[:-1])/values[1:]).tolist())
    (out/"orders.json").write_text(json.dumps(orders, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    run(parser.parse_args().output)
