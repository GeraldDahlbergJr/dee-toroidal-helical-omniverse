"""Exact certificate for one norm relaxation; no historical code is imported."""
import json
from pathlib import Path
import sympy as s

m = 6
x = s.Matrix(range(m))
one = s.ones(m, 1)
w = s.symbols('w0:6')
q = s.symbols('q01 q02 q03 q04 q05 q12 q13 q14 q15 q23 q24 q25 q34 q35 q45')
A = s.zeros(m)
A[0, 0] = -s.Rational(1, 2)
for z, (i, j) in zip(q, [(i,j) for i in range(m) for j in range(i+1,m)]):
    A[i,j] = z
    A[j,i] = -z
b = s.zeros(m,1)
b[5] = s.Rational(1,2)  # Q_56, fixed by centered row 6 and unit interior norm
assert A + A.T == s.diag(-1,0,0,0,0,0)
c = A*one+b
f = A*x+6*b-s.Matrix(w)
variables = list(q)+list(w)
# Only linear SBP/constant/affine constraints are eliminated.
point = next(iter(s.linsolve(list(c)+list(f), variables)))
sub = dict(zip(variables, point))
reduce = lambda z: s.factor(z.subs(sub, simultaneous=True))
assert all(reduce(z)==0 for z in list(c)+list(f))
t = A*x.applyfunc(lambda z:z*z)+36*b-2*s.diag(*w)*x
assert reduce(sum(w)-s.Rational(11,2)) == 0
assert reduce(sum(i*w[i] for i in range(m))-15) == 0
assert reduce(sum(t)-s.Rational(1,2)) == 0
# Coefficients in sum_{i=1}^6 i*w_i*R_i, including transition row 6.
# Q_65=-1/2 and x_6=6. The endpoint row has zero multiplier.
coeff = (x.T*A).T-6*b
assert all(reduce(coeff[j]+w[j])==0 for j in range(m))
# Polynomial system, valid iff rational moment system when every w_j>0.
product = s.prod(w)
P = []
for i in range(1,6):
    P.append(s.expand(sum(A[i,j]*t[j]*s.prod(w[k] for k in range(m) if k!=j) for j in range(m))))
P.append(-t[5]*s.prod(w[:5])/2)
# Factored algebra certifies the combination without expanding large reductions.
assert s.expand(sum(i*P[i-1] for i in range(1,7)) - sum(coeff[j]*t[j]*s.prod(w[k] for k in range(m) if k!=j) for j in range(m))) == 0
out = {
 'base':'fab11b7ed4defc34687d88ab7fdb277ee91998a7',
 'arithmetic':'exact rational/symbolic', 'sympy':s.__version__,
 'nodes':16, 'boundary_weights':[str(z) for z in w],
 'Q_inner_block':[[str(z) for z in A.row(i)] for i in range(m)],
 'constant_equations':[str(z) for z in c],
 'affine_equations':[str(z) for z in f],
 'quadratic_error_numerators':[str(s.expand(z)) for z in t],
 'all_affected_moment_equations':'P_i=sum_j Q_ij*t_j*prod_{k!=j}w_k=0, i=1,...,5; P_6=-t_5*prod_{k=0}^4 w_k/2=0; outer rows reflected',
 'forced_weight_sum':'11/2', 'forced_first_weight_moment':'15',
 'forced_weighted_error_sum':'1/2',
 'polynomial_contradiction':'sum_{i=1}^6 i*P_i = -prod_{j=0}^5 w_j/2; positivity makes RHS nonzero',
 'positive_weighted_norm':'INFEASIBLE', 'boundary_moments':'FAIL',
 'Dirichlet_strain_kernel':'UNRESOLVED', 'adjoint_identity':'UNRESOLVED',
 'stop_reason':'Exact necessary moment subsystem inconsistent; no member selected or checked further',
 'physical_solves':0, 'grid_sweeps':0, 'historical_certificates_rerun':False
}
Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(out['polynomial_contradiction'])
