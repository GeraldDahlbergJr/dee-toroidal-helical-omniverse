"""Compare completed historical and streamed-grid residual/field diagnostics."""
import json
import sys
from pathlib import Path
import numpy as np

historical,streamed=[json.loads(Path(path).read_text()) for path in sys.argv[1:]]
rtol=2e-5
atol=5e-13
comparisons=[]
for old,new in zip(historical['rows'][:2],streamed['rows'][:2]):
    assert old['n']==new['n']
    checks={}
    for q in (6,8):
        key=f'evaluator_q{q}'
        for metric in ('field_rms','field_max','residual_rms'):
            a,b=old[key][metric],new[key][metric]
            metrics=a.keys() if isinstance(a,dict) else (None,)
            for region in metrics:
                x=np.asarray(a if region is None else a[region])
                y=np.asarray(b if region is None else b[region])
                np.testing.assert_allclose(x,y,rtol=rtol,atol=atol)
                label=f'{key}/{metric}' if region is None else f'{key}/{metric}/{region}'
                checks[label]=float(np.max(abs(x-y)))
    for trace in ('radial_max','seam_max'):
        assert old['traces'][trace]==new['traces'][trace]
    assert new['linear_relative_max']<=1e-10
    comparisons.append({'n':old['n'],'max_abs':checks})
assert len(comparisons)==2
print(json.dumps({'rtol':rtol,'atol':atol,'comparisons':comparisons},indent=2))
