"""Independent exact Machin-pi and all denominator checks."""
from fractions import Fraction as F
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parents[3]

def atan_bounds(q,n=50):
    s=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n)),F(0))
    t=F((-1)**n,(2*n+1)*q**(2*n+1))
    return min(s,s+t),max(s,s+t)
a,b=atan_bounds(5);c,d=atan_bounds(239)
plo,phi=16*a-4*d,16*b-4*c
claimed=F('3.14159265358979323846264338327950288419716939937510')
assert claimed<=plo<=phi<=claimed+F(1,10**50)
rows=json.loads((ROOT/'research/run43_final_certification/error_budget.json').read_text())['scores']
for row in rows:
    n=row['n']
    with np.load(ROOT/f'research/run43_bspline_completion/coefficients_{n}.npz') as z:
        coeff=z['coefficients'][0]
    low=1+min(F(0),F(float(coeff.min())));high=1+max(F(0),F(float(coeff.max())))
    region=row['region']
    radial={'inner':F(3,800),'bulk':F(3,200),'outer':F(9,800),'first_inner':F(1,20)+F(1,5*n),'first_outer':F(1,4)-F(1,5*n)}[region]
    lo,hi=map(F,row['continuous_weight_denominator_exact_rational_enclosure'])
    assert lo<=4*plo**2*radial*low**6 and hi>=4*phi**2*radial*high**6
print(json.dumps({'pi_bracket_proven_by':'Machin identity and 50-term alternating rational arctangent series','denominator_rows_checked':len(rows),'domain':'exact decimal protocol chart','denominator_proof':'Nonnegative spline partition gives psi bounds; integrate positive J=r(1+r cos t) against these bounds; angular cos integral vanishes. Surface normalization uses protocol chart weight J rather than geometric area.','enclosure_arithmetic':'exact rational','status':'PASS for denominator only'},indent=2))
