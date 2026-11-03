"""Independent process: preserved frozen audit plus direct 4D source checks.

Does not change fields, source, audit, or solver output. Writes only the report
named on the command line. The frozen audit's endpoint handling is unchanged.
"""
import hashlib, importlib.util, json, math, pathlib, sys
import numpy as np

source = pathlib.Path(sys.argv[1]).resolve()
directory = pathlib.Path(sys.argv[2]).resolve()
destination = pathlib.Path(sys.argv[3]).resolve()
sys.path.insert(0, str(source / 'reproduce'))
from diagnose_toroidal_boundary_convergence import Audit

audit_path = source / 'reproduce/diagnose_toroidal_boundary_convergence.py'
raw = audit_path.read_bytes()
blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
assert blob == '5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
rows=[]
for n in (8,16,32,64):
    path=directory/f'fields_{n}.npz'
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == json.loads((directory/f'result_{n}.json').read_text())['sha256']
    with np.load(path) as data:
        psi=data['psi']; w=data['W_cartesian']
    assert psi.shape==(n+1,n,n) and w.shape==(3,n+1,n,n)
    assert np.isfinite(psi).all() and np.isfinite(w).all() and psi.min()>0
    assert np.array_equal(psi[[0,-1]],np.ones_like(psi[[0,-1]]))
    assert np.count_nonzero(w[:,[0,-1]])==0
    audit=Audit(psi.shape)
    # Reconstruct T_{mu nu} directly from the original action at actual solved
    # metric samples. No production projection or tensor-check helper is used.
    kinetic=np.array([[1,0,0],[0,1+.3*.945,.1*.945],[0,.1*.945,1+.2*.945]])
    pi=np.array([0.,.12,-.04]); worst_e=0.;worst_s=0.
    for i,j,k in [(0,0,0),(n//2,n//4,n//4),(n,n//2,n//2),(1,n//3,n//5)]:
        r=.05+.2*i/n;t=2*np.pi*j/n;p=2*np.pi*k/n;q=1+r*np.cos(t)
        et=np.array([-np.sin(t)*np.cos(p),-np.sin(t)*np.sin(p),np.cos(t)])
        ep=np.array([-np.sin(p),np.cos(p),0.])
        direction=et/r-3*ep/q
        derivative=np.zeros((3,4)); derivative[:,0]=pi
        derivative[1,1:]=.25*direction;derivative[2,1:]=.20*direction
        metric=np.diag([-1.,*[psi[i,j,k]**4]*3])
        contraction=np.einsum('ab,am,mn,bn',kinetic,derivative,np.linalg.inv(metric),derivative)
        lagrangian=-.5*contraction-(.945-1)**2
        stress=np.einsum('ab,am,bn->mn',kinetic,derivative,derivative)+metric*lagrangian
        expected_e=.5*audit.temporal+.5*audit.spatial[i,j,k]*psi[i,j,k]**-4+(.945-1)**2
        worst_e=max(worst_e,abs(stress[0,0]-expected_e))
        worst_s=max(worst_s,float(np.max(abs(-stress[1:,0]-audit.S[:,i,j,k]))))
    regions=audit.summarize(psi,w)
    rows.append(dict(n=n,sha256=digest,psi_min=float(psi.min()),psi_max=float(psi.max()),
                     direct_action_source_max_abs=dict(energy=worst_e,momentum=worst_s),regions=regions))
    print(json.dumps(dict(completed=n,sha256=digest,source_energy_difference=worst_e,source_momentum_difference=worst_s)),flush=True)
checks=[]
for band in ('inner_boundary_band','outer_boundary_band'):
    for metric in ('independent_flux_H','independent_tensor_M_order2'):
        for a,b in zip(rows,rows[1:]):
            order=float(np.log2(a['regions'][band][metric]['volume_weighted_rms']/b['regions'][band][metric]['volume_weighted_rms']))
            checks.append(dict(band=band,metric=metric,pair=[a['n'],b['n']],order=order,passed=bool(np.isfinite(order) and order>=1.5)))
assert len(checks)==12
report=dict(audit_git_blob=blob,passed=sum(c['passed'] for c in checks),total=12,
            status='PASS' if all(c['passed'] for c in checks) else 'FAIL',checks=checks,rows=rows,
            scope='Fresh array validation of frozen initial-constraint gate; no full matter evolution, absolute continuum error bound, or physical survival claim.')
destination.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
