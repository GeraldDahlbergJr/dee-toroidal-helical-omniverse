"""Read archived fields, reproduce all 60 q8 scores; never call solve()."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import sys,json,hashlib,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]/'run43_bspline_completion'
sys.path.insert(0,str(ROOT))
from evaluator import evaluate
start=time.monotonic()
archive=json.loads((ROOT/'results.json').read_text())
checks={}
for line in (ROOT/'SHA256SUMS').read_text().splitlines():
    expected,name=line.split(); actual=hashlib.sha256((ROOT/name).read_bytes()).hexdigest(); checks[name]={'expected':expected,'actual':actual,'match':expected==actual}
assert all(x['match'] for x in checks.values())
out={'kind':'empirical_reproduction_not_certificate','hashes':checks,'rows':[]}
for row in archive['rows']:
    n=row['n']; z=np.load(ROOT/f'coefficients_{n}.npz',allow_pickle=False); c=z['coefficients']
    assert int(z['n'])==n and c.shape==(4,n+2,n,n) and np.isfinite(c).all()
    k=int(os.environ.get("CERT_RULE","8")); now=evaluate(c,n,k)
    scores=[]
    for region,values in now['residual_rms'].items():
        for i,value in enumerate(values):
            old=row[f'evaluator_q{k}']['residual_rms'][region][i]
            scores.append({'region':region,'component':['H','Mx','My','Mz'][i],'archived':old,'reproduced':value,'absolute_difference':abs(value-old),'relative_difference':abs(value-old)/old})
    out['rows'].append({'n':n,'coefficient_shape':list(c.shape),'coefficient_min':c.min(axis=(1,2,3)).tolist(),'coefficient_max':c.max(axis=(1,2,3)).tolist(),'scores':scores,'elapsed_s':time.monotonic()-start})
    Path(__file__).with_name(f'reproduction_q{k}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(n, max(s['relative_difference'] for s in scores),flush=True)
