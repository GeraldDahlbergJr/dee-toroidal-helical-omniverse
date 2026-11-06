"""Read-only independent action/source and archived boundary checks; no PDE solve."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'reproduce'))
from dee_stress_energy_3p1 import projections
rng=np.random.default_rng(4303)
errors=np.zeros(4)
for _ in range(512):
    rho=rng.uniform(.2,1.5); p=rng.normal(size=3); D=rng.normal(size=(3,3))
    z=rng.normal(size=(3,3)); gamma=z.T@z+np.eye(3)
    alpha=rng.uniform(.5,2); beta=rng.normal(size=3)
    g=np.empty((4,4));g[0,0]=-alpha**2+beta@gamma@beta
    g[0,1:]=gamma@beta;g[1:,0]=g[0,1:];g[1:,1:]=gamma
    n=np.r_[1.,-beta]/alpha
    deriv=np.column_stack((alpha*p+D@beta,D))
    K=np.array([[1.,0,0],[0,1+.3*rho**2,.1*rho**2],[0,.1*rho**2,1+.2*rho**2]])
    lm=-.5*np.einsum('ab,am,mn,bn',K,deriv,np.linalg.inv(g),deriv)-(rho**2-1)**2
    T=deriv.T@K@deriv+g*lm
    direct=(n@T@n,-T[1:,:]@n,T[1:,1:],np.sum(np.linalg.inv(gamma)*T[1:,1:]))
    prod=projections(rho,p,D,np.linalg.inv(gamma))
    errors=np.maximum(errors,[np.max(np.abs(a-b)) for a,b in zip(direct,prod)])
P=np.array([.12,-.04]);d=np.array([.25,.20]);K=np.array([[1+.3*.945,.1*.945],[.1*.945,1+.2*.945]])
rows=[]
for ngrid in (8,16,32,64):
    path=ROOT/f'research/results/coordinated-replay-20261103-1831/baseline/fields_{ngrid}.npz'
    with np.load(path) as data: u=data['psi'];w=data['W_cartesian']
    rows.append(dict(n=ngrid,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),psi_min=float(u.min()),psi_max=float(u.max()),boundary_psi_max=float(np.max(abs(u[[0,-1]]-1))),boundary_W_max=float(np.max(abs(w[:,[0,-1]]))),finite=bool(np.isfinite(u).all() and np.isfinite(w).all())))
out=dict(seed=4303,samples=512,projection_max_abs=dict(zip(['E','S_i','S_ij','S_trace'],map(float,errors))),source_coefficients=dict(A=float(.5*P@K@P+(.945-1)**2),B_prefactor=float(.5*d@K@d),S_prefactor=float(-P@K@d),kinetic_eigenvalues=np.linalg.eigvalsh(K).tolist()),archived_boundaries=rows,scope='Algebra and stored nodal boundary checks only; no continuum residual or PDE solve.')
Path(__file__).with_name('equation_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
