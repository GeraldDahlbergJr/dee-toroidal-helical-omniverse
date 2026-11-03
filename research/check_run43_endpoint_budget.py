"""Polynomial checks for the endpoint derivation; no frozen data writes."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reproduce'))
from diagnose_toroidal_boundary_convergence import radial_matrix
from investigate_toroidal_inner_operator import endpoint_cubic_matrix

def main():
    checks=[]
    for n in (8,16,32,64):
        h=.2/n;x=np.arange(n+1)*h;d=radial_matrix(n+1,h,2)
        m=endpoint_cubic_matrix(n+1,h,True)
        for degree in (0,1,2,3,4):
            f=(x-x[1])**degree
            truth=2. if degree==2 else 0.
            expected=1.5*h if degree==3 else 5*h*h if degree==4 else 0.
            # Cubic: h/4*f''' = 1.5h. Quartic: 5h^2/24*f'''' = 5h^2.
            np.testing.assert_allclose((d@d@f)[1]-truth,expected,atol=1e-11,rtol=1e-10)
            exact=h/4*np.diff(f[:4],n=3)[0]/h**3
            np.testing.assert_allclose(((d@d-m@m)@f)[1],exact,atol=1e-11,rtol=1e-10)
            checks.append(dict(n=n,degree=degree,passed=True))
    inventory=json.loads(Path('/workspace/run43-verification/inventory.json').read_text())
    for e in inventory:
        assert hashlib.sha256((ROOT/e['path']).read_bytes()).hexdigest()==e['sha256']
    out=dict(polynomial_checks=checks,original_archive_hashes_unchanged=len(inventory),status='PASS')
    (ROOT/'research/results/endpoint-error-budget-20261103/checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('20 polynomial checks pass; all 105 original archive hashes unchanged')
if __name__=='__main__':main()
