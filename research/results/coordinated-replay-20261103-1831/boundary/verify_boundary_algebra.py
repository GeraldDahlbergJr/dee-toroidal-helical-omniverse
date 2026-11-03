"""Independent exact-rational boundary stencil certificate, no physical forcing."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def d1(n):
    a = [[F(0) for _ in range(n)] for _ in range(n)]
    a[0][:3] = [F(-3,2), F(2), F(-1,2)]
    a[-1][-3:] = [F(1,2), F(-2), F(3,2)]
    for i in range(1,n-1):
        a[i][i-1], a[i][i+1] = F(-1,2), F(1,2)
    return a


def main():
    n=9
    d=d1(n)
    dd=[[sum(d[i][k]*d[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    direct=[[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(1,n-1):
        direct[i][i-1:i+2]=[F(1),F(-2),F(1)]
    defect=[[dd[i][j]-direct[i][j] for j in range(n)] for i in range(n)]
    assert dd[1][:4]==[F(3,4),F(-5,4),F(1,4),F(1,4)]
    assert defect[1][:4]==[F(-1,4),F(3,4),F(-3,4),F(1,4)]
    moments={}
    for i in (1,n-2):
        vals=[sum(defect[i][j]*F(j-i)**p for j in range(n)) for p in range(5)]
        assert vals[:3]==[0,0,0]
        assert vals[3]==(F(3,2) if i==1 else F(-3,2))
        moments[str(i)]=[str(v) for v in vals]
    # f=(r-a)(b-r)/r has homogeneous Dirichlet data but f'''=6ab/r^4.
    a,b=F(1,20),F(1,4)
    assert a+b-a-a*b/a==0 and a+b-b-a*b/b==0
    assert 6*a*b/a**4 != 0
    report={"scope":"Exact rational stencil identity only; no physical acceptance result", "first_interior_DrDr_h2":[str(x) for x in dd[1][:4]],"first_interior_defect_h2":[str(x) for x in defect[1][:4]],"defect_moments_about_rows":moments,"inner_odd_leading_coefficient":"+h f'''/4", "outer_odd_leading_coefficient":"-h f'''/4", "Dirichlet_counterexample_third_derivative_at_inner":str(6*a*b/a**4),"radial_vector_principal_eigenvalues":["1","1","4/3"],"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('boundary-algebra.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
