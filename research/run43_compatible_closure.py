"""Original-source direct Hessian with consistent coordinate KO dissipation.

No coefficient is inferred from solved fields or acceptance orders.  The 1/4
coefficient is the identity D2 - D0 D0 = -h^2 D4/4 on a uniform interior.
Radial D4 uses degree-five polynomial closures, so no extra boundary datum
or modified matter source is introduced. All mixed Hessians are retained.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator
import run43_fourth_vector as fourth
import run43_consistent_momentum as baseline


def fourth_difference(n, h, periodic=False):
    count = 5 if periodic else 6
    a = sp.lil_matrix((n, n))
    for i in range(n):
        start = min(max(i-2, 0), n-count)
        js = np.arange(i-2, i+3) if periodic else np.arange(start, start+count)
        offsets = (js-i).astype(float)
        rhs = np.zeros(count); rhs[4] = 24
        weights = np.linalg.solve(np.array([offsets**k for k in range(count)]), rhs)/h**4
        for j, w in zip(js, weights):
            a[i, j % n] += w
    return a.tocsr()


class CompatibleGrid(fourth.FourthVectorGrid):
    def __init__(self, n, stabilization=.25):
        super().__init__(n, stabilization)
        eye = [sp.eye(k, format='csr') for k in self.shape]
        h = [.2/n, 2*np.pi/n, 2*np.pi/n]
        self.ko = []
        for axis in range(3):
            ops = eye.copy()
            ops[axis] = fourth_difference(self.shape[axis], h[axis], axis != 0)
            self.ko.append(sp.kron(sp.kron(ops[0], ops[1], format='csr'), ops[2], format='csr'))
        original = self.vector
        # Physical mesh spacing is a_a*h_a: (a_a*h_a)^2/a_a^4
        # gives h_a^2/a_a^2, rather than h_a^2/a_a^4.
        scale = [np.ones(self.size), self.r**-2, self.q**-2]
        def apply(x):
            full = np.zeros((3, self.size)); full[:, self.idx] = x.reshape(3, -1)
            result = original @ x
            correction = np.zeros_like(full)
            for axis in range(3):
                v = (self.ko[axis] @ full.T).T
                e = self.basis[axis]
                correction -= h[axis]**2/4 * scale[axis] * (v + e*np.sum(e*v, axis=0)/3)
            return result + correction[:, self.idx].ravel()
        self.vector = LinearOperator(original.shape, matvec=apply, dtype=float)


if __name__ == '__main__':
    solve = CompatibleGrid.solve
    def snapshot(self):
        out = baseline.output_dir()
        (out/'executed_compatible_closure.py').write_bytes(Path(__file__).read_bytes())
        (out/'executed_fourth_vector.py').write_bytes(Path(fourth.__file__).read_bytes())
        (out/'correction.json').write_text(json.dumps(dict(
            closure='direct fourth-order Cartesian Hessian; six-point radial polynomial D4 closure',
            dissipation='-sum_a h_a^2/4 scale_a^-2 (I+e_a e_a^T/3) D_aaaa W',
            coefficient_origin='centered D2-D0 squared identity; no fit',
            mixed='all direct Hessian mixed terms retained',
            source='original unchanged', frozen_audit='unchanged',
            extra_boundary_data=False), indent=2)+'\n')
        return solve(self)
    CompatibleGrid.solve = snapshot
    baseline.TensorGrid = CompatibleGrid
    baseline.main()
