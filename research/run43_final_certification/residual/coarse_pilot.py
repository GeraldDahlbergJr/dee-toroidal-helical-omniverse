"""Exact-rational, deliberately coarse analytic bound; see METHOD.md proof."""
from fractions import Fraction as Q
from pathlib import Path
import json,numpy as np
root=Path(__file__).resolve().parents[2]/'run43_bspline_completion'
n=32;c=np.load(root/'coefficients_32.npz')['coefficients']
assert np.isfinite(c).all() and np.max(np.abs(c))<=1 and np.min(c[0])>=0
# Verify every positive knot interval against conservative rational spacing.
A=float(.05);B=float(.25)
knots=[np.r_[[A]*5,np.linspace(A,B,n+1)[1:-1],[B]*5],np.arange(-4,n+5)*(2*np.pi/n)]
for k in knots:
 for a,b in zip(k[:-1],k[1:]):
  if b>a: assert Q(float(b))-Q(float(a))>=Q(1,10*n)
h=Q(1,10*n);D1=8/h;D2=64/h**2
H=4356*D2+99000*D1
Hex=4356*10000+99000*100
U=4*H+4*Hex+400000
a=json.loads((root/'results.json').read_text()); E=Q(a['rows'][2]['evaluator_q8']['residual_rms']['first_outer'][2])
out={'scope':'fixed saved field continuum RMS only; no exact-assembly solution claim','n':32,'region':'first_outer','component':'My','true_RMS_enclosure_rational':['0',str(U)],'absolute_difference_from_archived_score_bound_rational':str(max(E,U-E)),'absolute_difference_bound':float(max(E,U-E)),'allowed_total_error':float(E/10),'meets_10_percent':False,'method':'analytic pointwise range, positive weighted RMS; Fraction exact arithmetic','quadrature_truncation': 'bypassed by continuum pointwise range; no useful sharp truncation bound','solution_uncertainty':None}
Path(__file__).with_name('coarse_pilot.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
