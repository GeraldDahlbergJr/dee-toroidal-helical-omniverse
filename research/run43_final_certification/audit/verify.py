"""Independent exact-rational audit; no PDE solves and no empirical certification."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
SRC=ROOT/'research/run43_bspline_completion'
checks=[]
for directory in (SRC, ROOT/'research/run43_final_gate'):
    for line in (directory/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split(maxsplit=1); path=directory/name
        actual=hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
        checks.append(dict(path=str(path.relative_to(ROOT)),expected=expected,actual=actual,status='MATCH' if actual==expected else 'MISSING' if actual is None else 'MISMATCH'))
rows=json.loads((SRC/'results.json').read_text())['rows']
condition=json.loads((SRC/'conditioning.json').read_text())
coeff=[]
for row,archived in zip(rows,condition,strict=True):
    n=row['n']
    with np.load(SRC/f'coefficients_{n}.npz') as data:
        c=data['coefficients']; assert c.shape==(4,n+2,n,n) and c.dtype==np.float64 and np.isfinite(c).all()
        assert int(data['n'])==n
    # Recompute derivative bounds by EXACT subtraction of binary64 values.
    w=c[1:]; padded=np.pad(w,((0,0),(1,1),(0,0),(0,0)))
    diffs=[]
    for axis,values,periodic in ((1,padded,False),(2,w,True),(3,w,True)):
        moved=np.moveaxis(values,axis,0)
        pairs=zip(moved,np.roll(moved,1,axis=0)) if periodic else zip(moved[1:],moved[:-1])
        diffs.append(max(abs(F(float(a))-F(float(b))) for aa,bb in pairs for a,b in zip(aa.flat,bb.flat)))
    lo=1+min(F(0),F(float(c[0].min())));hi=1+max(F(0),F(float(c[0].max())))
    strain=6*(4*n/F('.20')*diffs[0]+4*n/(2*F('3.14159')*F('.05'))*diffs[1]+4*n/(2*F('3.14159')*F('.75'))*diffs[2])
    alpha=1-F('.034')*2*F('3.14160')*F('.001')*(5*F('.0127638')*hi**4+F('.068614375')*(400+9/F('.75')**2))
    v=F('.969');cross=F('.185')*strain/(2*lo**7);other=F('.034')*48*F('3.14160')*F('.001')*hi**5*F('.62');s=cross*other/(alpha*v)
    scalar=(1+cross/v)/(alpha*(1-s));vector=(1+other*scalar)/v
    assert F(archived['psi_lower'])<=lo and F(archived['psi_upper'])>=hi
    assert F(archived['strain_frobenius_upper'])>=strain and F(archived['scalar_gradient_coercivity_lower'])<=alpha
    assert F(archived['scalar_cross_gradient_norm_upper'])>=cross and F(archived['vector_cross_gradient_norm_upper'])>=other and F(archived['schur_gradient_norm_upper'])>=s
    # Separately reproduce final-gate bound from its outward archived inputs.
    a,vv,cc,dd,ss=[F(archived[k]) for k in ('scalar_gradient_coercivity_lower','vector_gradient_coercivity_lower','scalar_cross_gradient_norm_upper','vector_cross_gradient_norm_upper','schur_gradient_norm_upper')]
    independent_schur=cc*dd/(a*vv)
    assert independent_schur<1
    bs=(1+cc/vv)/(a*(1-independent_schur));bv=(1+dd*bs)/vv
    assert max(bs,bv)<F(133,100)
    coeff.append(dict(n=n,shape=list(c.shape),finite=True,exact_coefficient_difference_check=True,inverse_bound_exact=str(max(bs,bv)),inverse_bound_display=float(max(bs,bv))))
scores={(r['n'],region,comp):F(value) for r in rows for region,values in r['evaluator_q8']['residual_rms'].items() for comp,value in zip(('H','Mx','My','Mz'),values,strict=True)}
assert len(scores)==60
pairs=[]
for (n,region,comp),coarse in scores.items():
    if n==32:continue
    ratio=9*coarse/(11*scores[2*n,region,comp]);passed=ratio>0 and ratio*ratio>=8
    pairs.append(dict(n_coarse=n,n_fine=2*n,region=region,component=comp,conditional_ratio_exact=str(ratio),conditional_order_display=math.log2(float(ratio)),conditional_pass=passed,unconditional_pass=None))
assert len(pairs)==40 and all(x['conditional_pass'] for x in pairs)
pilot=scores[32,'first_outer','My']
out=dict(checkpoint='c5580eea5b216052ba5e91c6725e0304a7852a99',hash_checks=checks,coefficient_checks=coeff,score_count=len(scores),comparisons=pairs,minimum_conditional_order=min(x['conditional_order_display'] for x in pairs),pilot=dict(score=float(pilot),absolute_allowance_exact=str(pilot/10),absolute_allowance_display=float(pilot/10)),verdict='NOT YET ESTABLISHED',reason='No complete rigorous quadrature, floating-point and exact-assembly solution uncertainty enclosure supplied.')
json.dump(out,sys.stdout,indent=2);print()
