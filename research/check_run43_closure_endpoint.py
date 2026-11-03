"""Validate the endpoint discrepancy on a newly solved correction's fields."""
import sys,json
from pathlib import Path
import numpy as np
from validate_original_physics import Independent


def main():
    directory=Path(sys.argv[1]);target=directory/'endpoint-discrepancy.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    original=json.loads((directory/'original-physics-audit.json').read_text())
    for k,n in enumerate((8,16,32,64)):
        with np.load(directory/f'fields_{n}.npz') as z:u,w=z['psi'],z['W_cartesian']
        ref=Independent(n);h=.2/n
        delta3=np.diff(w[:,:4],n=3,axis=1)[:,0]/h**3;er=ref.er[:,1]
        endpoint=u[1]**-10*h/4*(delta3+er*np.sum(er*delta3,axis=0)/3)
        weight=np.broadcast_to(ref.r*ref.q,u.shape)*u**6
        mask=np.broadcast_to((ref.r>.05+1e-12)&(ref.r<=.10+1e-12),u.shape)
        norm=float(np.sqrt(np.sum(np.sum(endpoint**2,axis=0)*weight[1])/np.sum(weight[mask])))
        result=json.loads((directory/f'result_{n}.json').read_text())
        rows.append(dict(n=n,frozen_inner_M=result['audit']['inner_boundary_band']['independent_tensor_M_order2']['volume_weighted_rms'],exact_endpoint_difference_rms=norm,supplementary_inner_M=original['rows'][k]['regions']['inner']['M']['volume_weighted_rms']))
    target.write_text(json.dumps(dict(scope='Exact frozen-minus-matched endpoint difference on new fields only; not additive scalar budget or replacement gate.',rows=rows),indent=2)+'\n')
    print(target)


if __name__=='__main__':main()
