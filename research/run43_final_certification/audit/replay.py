"""Independent replay invocation of archived strong evaluator, no solves.
This checks reproducibility, not arithmetic enclosure or quadrature truncation.
"""
from pathlib import Path
import sys,json,time,numpy as np
ROOT=Path(__file__).resolve().parents[3]
SRC=ROOT/'research/run43_bspline_completion'
sys.path.insert(0,str(SRC))
from evaluator import evaluate
rows=json.loads((SRC/'results.json').read_text())['rows'];out=[]
for row in rows:
    n=row['n'];start=time.monotonic()
    with np.load(SRC/f'coefficients_{n}.npz') as z:c=z['coefficients']
    result=evaluate(c,n,8)['residual_rms']
    for region,values in result.items():
        for comp,actual,old in zip(('H','Mx','My','Mz'),values,row['evaluator_q8']['residual_rms'][region],strict=True):
            out.append(dict(n=n,region=region,component=comp,actual=actual,archived=old,absolute_difference=abs(actual-old),relative_difference=abs(actual-old)/old))
    print(f'n={n} replay seconds={time.monotonic()-start:.2f}',file=sys.stderr,flush=True)
assert len(out)==60
print(json.dumps(dict(scores=out,maximum_relative_difference=max(r['relative_difference'] for r in out),certified_error_bound=False),indent=2))
