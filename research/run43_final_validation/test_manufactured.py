"""Qualification gate; exits 1 on failure, without attempting repair or a solve."""
import json
import numpy as np
from evaluate import OUT, radial

def main():
    rows=[]
    for n in (8,16,32,64):
        matrix=radial(n,9,2)
        # Standard gamma_m bound for a floating-point dot product of the
        # STORED matrix with exact ones. No physical accuracy tolerance chosen.
        eps=np.finfo(float).eps
        gamma=(n+1)*eps/(1-(n+1)*eps)
        bound=gamma*np.max(np.sum(abs(matrix),axis=1))
        defect=np.max(abs(matrix@np.ones(n+1)))
        rows.append({'n':n,'constant_second_derivative_defect':float(defect),'stored_matrix_dot_product_roundoff_bound':float(bound),'passed':bool(defect<=bound)})
    # This necessary polynomial consistency check exposes coefficient-construction
    # error beyond matrix-application roundoff; it is not an error bound on all
    # arithmetic in coefficient construction or on the physical residual.
    fixture=json.loads((OUT/'results.json').read_text())['fixture']
    report={'status':'FAILED' if not all(r['passed'] for r in rows) else 'NOT YET ESTABLISHED','necessary_consistency_check':rows,'manufactured_fixture':fixture,'consequence':'Evaluator not qualified. Baseline continuum residuals, uncertainty and orders are provisional and cannot decide physical validation.'}
    (OUT/'qualification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    raise SystemExit(1)
if __name__=='__main__':main()
