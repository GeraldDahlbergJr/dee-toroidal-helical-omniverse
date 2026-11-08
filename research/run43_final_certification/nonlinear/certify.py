"""Reproduce sampled conditioning; never upgrades sampled evidence to exact assembly."""
import importlib.util,json,math,sys
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3]
SRC=ROOT/'research/run43_bspline_completion'
spec=importlib.util.spec_from_file_location('conditioning',SRC/'conditioning.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
archived=json.loads((SRC/'conditioning.json').read_text())
rows=[]
for n,old in zip((8,16,32),archived,strict=True):
    fresh=m.certify(SRC/f'coefficients_{n}.npz')
    assert fresh==old
    a,v,c,d,s=[F(fresh[k]) for k in ('scalar_gradient_coercivity_lower','vector_gradient_coercivity_lower','scalar_cross_gradient_norm_upper','vector_cross_gradient_norm_upper','schur_gradient_norm_upper')]
    # Compute Schur factor from outward inputs, independently of archived s.
    s_ind=c*d/(a*v)
    assert s_ind<1
    bs=(1+c/v)/(a*(1-s_ind));bv=(1+d*bs)/v
    cap=max(bs,bv)
    assert cap<F(133,100)
    rows.append({'n':n,'conditioning_reproduced':True,'independent_schur_exact':str(s_ind),'sampled_inverse_bound_exact':str(cap),'sampled_inverse_upper':math.nextafter(float(cap),math.inf),'exact_assembly_error_upper':None,'absolute_algebraic_dual_defect_upper':None,'jacobian_perturbation_upper':None,'nonlinear_remainder_constant_upper':None,'verified_solution_radius':None,'propagated_strong_residual_error_upper':None,'status':'NOT YET ESTABLISHED'})
print(json.dumps({'scope':'Mathematical sampled Gauss-6 Jacobian at exact archived binary64 coefficients; no new PDE solves','rows':rows},indent=2))
