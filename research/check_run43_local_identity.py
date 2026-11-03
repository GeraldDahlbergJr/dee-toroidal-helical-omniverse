"""Separate-process residual decomposition for local-boundary solved arrays."""
import sys
import json
from pathlib import Path
import numpy as np
from run43_local_boundary import LocalBoundaryGrid
from diagnose_toroidal_boundary_convergence import Audit


def main():
    source = Path(sys.argv[1])
    out = Path(sys.argv[2])
    out.mkdir(exist_ok=False)
    rows = []
    for n in (8, 16, 32, 64):
        with np.load(source/f'fields_{n}.npz') as z:
            u, w = z['psi'], z['W_cartesian']
        g = LocalBoundaryGrid(n)
        a = Audit(u.shape)
        flat = u.ravel()
        wf = w.reshape(3, -1)
        operator = (g.vector@wf[:, g.idx].ravel()).reshape(3, -1)
        rhs = 8*np.pi*.001*flat[g.idx]**6*g.S[:, g.idx]
        t = a.tensor(w, 2)
        tensor = np.array([sum(a.grad(t[i, j], 2)[j]
                              for j in range(3)) for i in range(3)])
        tensor = tensor.reshape(3, -1)[:, g.idx]
        algebraic = flat[g.idx]**-10*(operator-rhs)
        mismatch = flat[g.idx]**-10*(tensor-operator)
        frozen = a.momentum(u, w, 2).reshape(3, -1)[:, g.idx]
        closure = float(np.max(abs(frozen-algebraic-mismatch)))
        if closure > 1e-12:
            raise RuntimeError(f'Residual identity failed: {closure}')
        bands = {}
        for name, sl in [('first', slice(0, n*n)),
                         ('last', slice(-n*n, None))]:
            weight = (g.r*g.q*flat**6)[g.idx[sl]]
            def rms(v):
                return float(np.sqrt(np.sum(np.sum(v[:, sl]**2, axis=0)*weight)
                                     / np.sum(weight)))
            bands[name] = dict(algebraic=rms(algebraic),
                               mismatch=rms(mismatch), frozen=rms(frozen))
        rows.append(dict(n=n, identity_max_abs=closure, planes=bands))
        np.savez_compressed(out/f'planes_{n}.npz',
                            algebraic=algebraic[:, np.r_[0:n*n, -n*n:0]],
                            mismatch=mismatch[:, np.r_[0:n*n, -n*n:0]],
                            frozen=frozen[:, np.r_[0:n*n, -n*n:0]])
    report = dict(identity='M_frozen=psi^-10(A W-b)+psi^-10(T W-A W)',
                  source=str(source), rows=rows)
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
