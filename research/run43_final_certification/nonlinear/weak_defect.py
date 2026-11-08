"""No-solve binary64 diagnostics, NOT outward-rounded error certificates."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; SRC=ROOT/'research/run43_bspline_completion'
sys.path.insert(0,str(SRC))
import experiment as e
n=int(sys.argv[1]) if len(sys.argv)>1 else 32
c=np.load(SRC/f'coefficients_{n}.npz')['coefficients']
k=int(sys.argv[2]) if len(sys.argv)>2 else 6
start=time.monotonic(); space=e.Space(n,k)
# Reassemble weak defects at the stored coefficients, without a solve.
f=[e.forcing(b) for b in space.slabs()]
sq=[e.strain_squared(b,c[1:]) for b in space.slabs()]
rh,_=e.scalar_rhs(space,c[0],f,sq)
rm=e.vector_rhs(space,c[0],f)-space.vector(c[1:])
# K_g >= .0375 H on the common tensor space for physical chart bounds.
# pre() = H^-1/.15, so ||r||_(g*) <= sqrt(4 r.pre(r)) in exact arithmetic.
# This executable does NOT enclose assembly or spectral eigensolve roundoff.
values=[float(np.sqrt(max(0,4*np.vdot(r,space.pre(r))))) for r in [rh,*rm]]
result={'n':n,'quadrature':k,'seconds':time.monotonic()-start,'scalar_dual_upper_formula_binary64':values[0],'vector_dual_upper_formula_binary64':float(np.linalg.norm(values[1:])),'coefficient_defect_euclidean_norms':[float(np.linalg.norm(rh)),float(np.linalg.norm(rm))],'rigorous':False,'reason':'Floating-point assembly, symbolic-source evaluation, and eigensolve/preconditioner not enclosed; exact-minus-Gauss assembly absent.'}
if '--compare' in sys.argv:
    other=e.Space(n,8 if k==6 else 6)
    ff=[e.forcing(b) for b in other.slabs()]
    ss=[e.strain_squared(b,c[1:]) for b in other.slabs()]
    hh,_=e.scalar_rhs(other,c[0],ff,ss)
    mm=e.vector_rhs(other,c[0],ff)-other.vector(c[1:])
    dd=[float(np.sqrt(max(0,4*np.vdot(r,space.pre(r))))) for r in [hh-rh,*(mm-rm)]]
    result['q8_minus_q6_dual_formula_binary64']={'scalar':dd[0],'vector':float(np.linalg.norm(dd[1:]))}
print(json.dumps(result,indent=2))
