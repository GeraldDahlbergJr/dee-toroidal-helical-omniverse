"""New row-local residual identity on existing original-source fourth fields.

No solve or acceptance audit is repeated. Separates algebraic solve residual
from operator mismatch at the row under investigation, including psi factors.
"""
import sys,json
from pathlib import Path
import numpy as np
from run43_fourth_vector import FourthVectorGrid
from diagnose_toroidal_boundary_convergence import Audit


def main():
    out=Path(sys.argv[1]);source=Path(sys.argv[2])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);rows=[]
    for n in (8,16,32,64):
        with np.load(source/f'fields_{n}.npz') as z:u,w=z['psi'],z['W_cartesian']
        g=FourthVectorGrid(n);a=Audit(u.shape);flat=u.ravel();wf=w.reshape(3,-1)
        H=(g.vector@wf[:,g.idx].ravel()).reshape(3,-1)
        b=8*np.pi*.001*flat[g.idx]**6*g.S[:,g.idx]
        t=a.tensor(w,2)
        T=np.array([sum(a.grad(t[i,j],2)[j] for j in range(3)) for i in range(3)]).reshape(3,-1)[:,g.idx]
        algebraic=flat[g.idx]**-10*(H-b)
        mismatch=flat[g.idx]**-10*(T-H)
        frozen=a.momentum(u,w,2).reshape(3,-1)[:,g.idx]
        closure=float(np.max(abs(frozen-algebraic-mismatch)))
        assert closure<1e-12
        def rownorm(v):
            v=v[:,:n*n];weight=(g.r*g.q*flat**6)[g.idx[:n*n]]
            return float(np.sqrt(np.sum(np.sum(v*v,axis=0)*weight)/np.sum(weight)))
        row=dict(n=n,first_row_radius=float(g.r[g.idx[0]]),algebraic_first_row_rms=rownorm(algebraic),mismatch_first_row_rms=rownorm(mismatch),frozen_first_row_rms=rownorm(frozen),identity_max_abs=closure)
        rows.append(row)
        np.savez_compressed(out/f'row_{n}.npz',solver_physical_residual=algebraic[:,:n*n],operator_mismatch=mismatch[:,:n*n],frozen_residual=frozen[:,:n*n])
    (out/'report.json').write_text(json.dumps(dict(source=str(source),identity='M_frozen = psi^-10 (H W-b) + psi^-10 (T W-H W)',scope='Original-source saved fields, row-local identity only; previous all-twelve audit remains frozen.',rows=rows),indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
