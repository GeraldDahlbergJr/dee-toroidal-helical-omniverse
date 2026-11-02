"""Recompute all 12 frozen checks from raw arrays in a separate process."""
import sys,json,hashlib,subprocess
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'reproduce'))
from diagnose_toroidal_boundary_convergence import Audit

def main():
    directory=Path(sys.argv[1]);audit=root/'reproduce/diagnose_toroidal_boundary_convergence.py'
    blob=subprocess.check_output(['git','hash-object',str(audit)],text=True).strip()
    assert blob=='5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
    rows=[]
    for n in (8,16,32,64):
        path=directory/f'fields_{n}.npz';raw=path.read_bytes()
        expected=json.loads((directory/f'result_{n}.json').read_text())['sha256']
        assert hashlib.sha256(raw).hexdigest()==expected
        with np.load(path) as z:u,w=z['psi'],z['W_cartesian']
        assert u.shape==(n+1,n,n) and w.shape==(3,n+1,n,n)
        assert np.isfinite(u).all() and np.isfinite(w).all() and u.min()>0
        assert np.array_equal(u[[0,-1]],np.ones_like(u[[0,-1]]))
        assert np.count_nonzero(w[:,[0,-1]])==0
        rows.append(dict(n=n,sha256=expected,audit=Audit(u.shape).summarize(u,w)))
    checks=[]
    for band in ('inner_boundary_band','outer_boundary_band'):
        for metric in ('independent_flux_H','independent_tensor_M_order2'):
            for a,b in zip(rows,rows[1:]):
                order=float(np.log2(a['audit'][band][metric]['volume_weighted_rms']/b['audit'][band][metric]['volume_weighted_rms']))
                checks.append(dict(band=band,metric=metric,pair=[a['n'],b['n']],order=order,passed=bool(np.isfinite(order) and order>=1.5)))
    report=dict(audit_git_blob=blob,passed=sum(c['passed'] for c in checks),total=12,status='PASS' if all(c['passed'] for c in checks) else 'FAIL',checks=checks,rows=rows,scope='Independent array re-evaluation of frozen numerical gate; not independent scientific replication or continuum proof.')
    (directory/'independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
    if report['status']!='PASS':sys.exit(1)
if __name__=='__main__':main()
