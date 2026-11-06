"""Verify primitive arithmetic bounds and evidence decisions without solves."""
from pathlib import Path
import sys,json
import numpy as np
O=Path(__file__).resolve().parent
sys.path.insert(0,str(O.parent/'run43_final_validation'))
from evaluate import radial
q=json.loads((O/'qualification.json').read_text());r=json.loads((O/'results.json').read_text())
assert q['status']=='PASS' and r['qualification_status'].startswith('PASS:')
assert [a['n'] for a in r['rows']]==[8,16,32,64]
for a,b,order in zip(q['rows'][1:],q['rows'][2:],q['orders']):
    ac=np.array(a['max_abs_error']);bc=np.array(b['max_abs_error']);ua=np.array(a['floating_error_bound']);ub=np.array(b['floating_error_bound'])
    assert np.all(ac>10*ua) and np.all(bc>10*ub)
    low=np.log2((ac-ua)/(bc+ub))
    assert np.all(low>=1.5) and np.allclose(low,order['lower_order'],rtol=1e-14,atol=0)
checks=[];eps=np.finfo(float).eps
for n in (8,16,32,64):
    gamma=(n+1)*eps/(1-(n+1)*eps)
    for d in (1,2):
        matrix=radial(n,9,d);bound=gamma*np.max(np.sum(abs(matrix),axis=1));error=np.max(abs(matrix@np.ones(n+1)))
        assert error<=bound
        checks.append(dict(n=n,derivative=d,test='constant',error=float(error),bound=float(bound),passed=True))
        for k in range(9):
            s=np.arange(n+1,dtype=float)/n;values=s**k;factor=1
            for j in range(d):factor*=k-j
            exact=factor*5**d*s**max(k-d,0) if k>=d else np.zeros(n+1)
            bounds=(eps/(1-eps)+gamma)*(abs(matrix)@abs(values))
            errors=abs(matrix@values-exact);assert np.all(errors<=bounds)
            checks.append(dict(n=n,derivative=d,test='polynomial',power=k,max_error=float(errors.max()),max_bound=float(bounds.max()),passed=True))
ids=json.loads((O/'continuum_identity_checks.json').read_text());assert ids['status']=='PASS' and len(ids['checks'])==12
pres=json.loads((O/'preservation.json').read_text());assert pres['all_original_tracked_blobs_unchanged'] and all(pres['original_input_hash_checks'].values())
assert all(v['matches_saved_result'] for v in pres['accepted_archives'].values())
for row in r['rows']:assert row['boundary']['psi_trace_max_error']==row['boundary']['W_trace_max_error']==0
assert any(not all(v['passes_p_1_5']) for v in r['orders'])
assert any(not all(z['within_10_percent']) for row in r['rows'] for z in row['regions'].values())
for p in O.glob('*.json'):json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
out={'status':'PASS','checks':'Qualification lower orders recomputed; 8/8 constant and 72/72 polynomial tests pass unchanged; 12 symbolic identities; archive identities and radial traces; failed baseline gate retained; JSON finite.','primitive_checks':checks,'publication_readiness':'NOT READY','new_DEE_solves':0}
(O/'verification.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(out['checks'])
