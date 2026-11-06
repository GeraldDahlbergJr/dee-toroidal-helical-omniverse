"""Read-only baseline endpoint budget; new purpose: isolate accepted baseline cancellation.

No solves and no acceptance substitution. Run from repository root with
OPENBLAS_NUM_THREADS=1 python research/priority_one_validation/agent1/diagnose_baseline.py
"""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'reproduce'),str(ROOT/'research')]
from diagnose_toroidal_boundary_convergence import Audit
from investigate_toroidal_inner_operator import endpoint_cubic_matrix

rows=[]
for n in (8,16,32,64):
    path=ROOT/f'research/results/coordinated-replay-20261103-1831/baseline/fields_{n}.npz'
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest==json.loads(path.with_name(f'result_{n}.json').read_text())['sha256']
    with np.load(path) as d:u,w=d['psi'],d['W_cartesian']
    a=Audit(u.shape);f=a.momentum(u,w,2);fourth=a.momentum(u,w,4)
    b=Audit(u.shape);b.drad[2]=endpoint_cubic_matrix(n+1,a.dr,True)
    rest=b.momentum(u,w,2);endpoint=f-rest
    weight=np.broadcast_to(a.r*a.q,u.shape)*u**6
    mask=np.broadcast_to((a.r>.05+1e-12)&(a.r<=.1+1e-12),u.shape)
    def dot(v,z):return float(np.sum((np.sum(v*z,axis=0)*weight)[mask])/np.sum(weight[mask]))
    def norm(v):return float(np.sqrt(max(0,dot(v,v))))
    d3=np.diff(w[:,:4],n=3,axis=1)[:,0]/a.dr**3;er=a.er[:,1]
    exact=np.zeros_like(w);exact[:,1]=u[1]**-10*a.dr/4*(d3+er*np.sum(er*d3,axis=0)/3)
    assert norm(endpoint-exact)<1e-12
    def components(v):return [float(np.sqrt(np.sum((v[c]**2*weight)[mask])/np.sum(weight[mask]))) for c in range(3)]
    row=dict(n=n,sha256=digest,frozen_rms=norm(f),fourth_rms=norm(fourth),endpoint_rms=norm(endpoint),matched_rms=norm(rest),endpoint_formula_relative_error=norm(endpoint-exact)/norm(endpoint),endpoint_remainder_cosine=dot(endpoint,rest)/norm(endpoint)/norm(rest),cancellation_fraction=-2*dot(endpoint,rest)/(dot(endpoint,endpoint)+dot(rest,rest)),frozen_components=components(f),fourth_components=components(fourth),difference_rms=norm(f-fourth))
    def d4(v,axis):return (np.roll(v,2,axis)-4*np.roll(v,1,axis)+6*v-4*np.roll(v,-1,axis)+np.roll(v,-2,axis))/(2*np.pi/n)**2
    stabilizer=np.array([.25*(d4(v,1)/a.r**2+d4(v,2)/a.q**2) for v in w])*u**-10
    row['solver_stabilizer_rms']=norm(stabilizer)
    row['frozen_minus_stabilizer_rms']=norm(f-stabilizer)
    rows.append(row)
    print(json.dumps(row),flush=True)
pairs=[]
for a,b in zip(rows,rows[1:]):
    pairs.append(dict(pair=[a['n'],b['n']],orders={k:float(np.log2(a[k]/b[k])) for k in ['frozen_rms','fourth_rms','endpoint_rms','matched_rms']},fourth_component_orders=[float(np.log2(x/y)) for x,y in zip(a['fourth_components'],b['fourth_components'])]))
out=dict(scope='Accepted baseline saved fields only. Closure subtraction and fourth derivatives are supplementary diagnostics, never replacement acceptance tests.',rows=rows,pairs=pairs)
Path(__file__).with_name('diagnostics.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(pairs,indent=2))
