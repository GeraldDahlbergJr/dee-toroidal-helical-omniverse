"""Exact rational error targets, NOT a quadrature certificate.
All arithmetic for displayed derivative ceilings uses Fraction; pi<=22/7.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
root=Path(__file__).resolve().parents[2]/'run43_bspline_completion'
a=json.loads((root/'results.json').read_text())
r=next(x for x in a['rows'] if x['n']==32)
e=F(r['evaluator_q8']['residual_rms']['first_outer'][2]); e6=F(r['evaluator_q6']['residual_rms']['first_outer'][2])
C=F(factorial(8)**4,17*factorial(16)**3)
# For first outer surface r=1/4-(1/5)/32=39/160, q>=3/4, u>=1.
# L=2*pi >=6, h=2*pi/32 <=11/56.
Dlower=F(39,160)*F(3,4)*36
Lupper=F(44,7); hupper=F(11,56)
# Allocate 0.10 E^2 Dlower to numerator truncation alone. The full
# 10% RMS requirement remains unproved even if this allocation is met.
Ntarget=F(1,10)*e*e*Dlower
ceiling=Ntarget/(Lupper**2*C*hupper**16)
ratio=abs(e-e6)/e
out={'status':'NOT YET ESTABLISHED','n':32,'region':'first_outer','component':'My','archived_q8':float(e),'absolute_total_allowance':float(e/10),'empirical_q6_q8_difference':float(abs(e-e6)),'empirical_relative_difference':float(ratio),'surface_D_lower_bound_rational':str(Dlower),'gauss8_remainder_constant_rational':str(C),'illustrative_numerator_truncation_budget':float(Ntarget),'sufficient_sum_of_16th_derivative_suprema_ceiling':float(ceiling),'ceiling_exact_rational':str(ceiling),'certified_derivative_suprema':None,'roundoff_bound':None,'propagated_solution_bound':None}
Path(__file__).with_name('pilot.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
