"""One-shot qualification gate. No baseline evaluation, repair loop or solve."""
import json
import numpy as np
from evaluate import OUT, radial, fixture

def main():
    rows=[]; polynomials=[]
    eps=np.finfo(float).eps
    for n in (8,16,32,64):
        gamma=(n+1)*eps/(1-(n+1)*eps)
        for derivative in (1,2):
            matrix=radial(n,9,derivative)
            bound=gamma*np.max(np.sum(abs(matrix),axis=1))
            defect=np.max(abs(matrix@np.ones(n+1)))
            rows.append(dict(n=n,derivative=derivative,defect=float(defect),
                             stored_matrix_dot_product_roundoff_bound=float(bound),passed=bool(defect<=bound)))
            # s=i/n is dyadic on these grids, as are all s**k, k<=8.
            # Exact derivative is k!/(k-d)! * 5**d * s**(k-d).
            s=np.arange(n+1,dtype=float)/n
            for k in range(9):
                values=s**k
                factor=1
                for j in range(derivative): factor*=k-j
                exact=factor*5**derivative*s**max(k-derivative,0) if k>=derivative else np.zeros(n+1)
                actual=matrix@values
                # Rounded exact rational coefficients plus stored dot product.
                bounds=(eps/(1-eps)+gamma)* (abs(matrix)@abs(values))
                error=abs(actual-exact)
                polynomials.append(dict(n=n,derivative=derivative,power=k,max_abs_error=float(error.max()),
                                        max_roundoff_bound=float(bounds.max()),passed=bool(np.all(error<=bounds))))
    manufactured=fixture()
    orders=[]
    for a,b in zip(manufactured['rows'][1:],manufactured['rows'][2:]):
        p=np.log2(np.array(a['max_abs_error'])/b['max_abs_error'])
        orders.append(dict(pair=[a['n'],b['n']],order=p.tolist(),passes_p_1_5=(p>=1.5).tolist()))
    passed=all(r['passed'] for r in rows+polynomials) and all(all(r['passes_p_1_5']) for r in orders)
    report=dict(status='PASS' if passed else 'FAIL',components=['H','Mx','My','Mz'],
                necessary_consistency_check=rows,polynomial_derivative_checks=polynomials,
                manufactured_fixture=manufactured,finest_pair_orders=orders,
                criteria='Constant application roundoff bound unchanged. Polynomial bound includes final coefficient rounding and dot-product rounding; exact dyadic test inputs. Manufactured p>=1.5 on 16->32 and 32->64; no fitted absolute tolerance or roundoff exemption.',
                consequence='Qualification failed: STOP. Do not evaluate baseline or initiate another repair.' if not passed else 'Qualified for subsequent baseline assessment.')
    (OUT/'repair_qualification.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('polynomial_derivative_checks','manufactured_fixture')},indent=2))
    raise SystemExit(0 if passed else 1)
if __name__=='__main__':main()
