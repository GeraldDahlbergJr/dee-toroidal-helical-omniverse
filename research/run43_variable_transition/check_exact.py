"""One variable-transition family; exact algebra only, no solver imports."""
import json
from pathlib import Path
import sympy as s

M = 8
H = s.diag(s.Rational(1,2), *([1]*(M-1)))
x = s.Matrix(range(M))
one = s.ones(M,1)
A = s.zeros(M)
A[0,0] = -s.Rational(1,2)
q = []
for i in range(M):
    for j in range(i+1,M):
        z = s.Symbol(f'q{i}{j}')
        q.append(z)
        A[i,j], A[j,i] = z, -z
b = s.zeros(M,1)
b[M-1] = s.Rational(1,2)
c = A*one+b
f = A*x+M*b-H*one
assert A+A.T == s.diag(-1,*([0]*(M-1)))
# Linear constant/affine elimination only; no nonlinear equation solving.
point = next(iter(s.linsolve(list(c)+list(f), q)))
sub = dict(zip(q,point))
A0 = A.subs(sub, simultaneous=True)
assert A0*one+b == s.zeros(M,1)
assert A0*x+M*b == H*one
e = H.inv()*(A0*x.applyfunc(lambda z:z*z)+M*M*b)-2*x
Rblock = H.inv()*A0*e
R = [s.factor(Rblock[i]) for i in range(1,M)] + [-s.factor(e[M-1])/2]
weighted_error = s.factor((one.T*H*e)[0])
weighted_residual = s.factor(sum(i*R[i-1] for i in range(1,M+1)))
assert s.simplify((x.T*A0-M*b.T)+one.T*H) == s.zeros(1,M)

# Independent exact abstract weighted-adjoint check, performed before reporting
# the moment inconsistency. These symbolic blocks represent arbitrary L0, Ht,
# and arbitrary boundary traction pairing; they do not model a physical grid.
ell = s.Matrix(2,2,s.symbols('l0:4'))
ht = s.diag(*s.symbols('t0:2', positive=True))
mv = s.diag(*s.symbols('v0:2', positive=True))
G = ell.T*ht*ell/2
AD = -mv.inv()*G
assert mv*AD == -G
assert G == G.T
u = s.Matrix(s.symbols('u0:2'))
assert s.expand((u.T*G*u)[0]-(ell*u).dot(ht*(ell*u))/2) == 0

# Final feasibility decision and mandatory stop: no kernel investigation follows.
assert weighted_error == s.Rational(1,2)
assert weighted_residual == -s.Rational(1,2)

free = sorted(A0.free_symbols,key=str)
out = {
 'base':'51bfca8feb707f0716ec3efc23a486389e80eb5d',
 'family':'six boundary rows plus two variable transition rows; arbitrary skew Q on nodes 0..7; Q78=1/2; centered from row 8',
 'arithmetic':'exact rational/symbolic', 'sympy':s.__version__,
 'largest_assembled_matrix':[8,8], 'physical_grids_assembled':0,
 'free_parameters_after_linear_constraints':[str(z) for z in free],
 'constant_equations':[str(z) for z in c],
 'affine_equations':[str(z) for z in f],
 'Q_block_after_linear_elimination':[[str(s.factor(A0[i,j])) for j in range(M)] for i in range(M)],
 'quadratic_first_derivative_errors':[str(s.factor(z)) for z in e],
 'all_quadratic_residuals_rows_1_through_8':[str(z) for z in R],
 'weighted_error_sum':str(weighted_error),
 'exact_contradiction':'sum_{i=1}^8 i*R_i = -1/2',
 'weighted_residual_sum':str(weighted_residual),
 'scalar_SBP':'PASS',
 'adjoint_identity':'PASS: formal full STF transpose/traction identity in report; abstract Gram algebra checked here',
 'positive_energy':'PASS: positive semidefinite Gram form; no strict positivity claimed',
 'boundary_moments':'INFEASIBLE', 'Dirichlet_strain_kernel':'UNRESOLVED: stopped at inconsistent moments',
 'coercivity':'UNRESOLVED', 'physical_solves':0, 'GMRES':0,
 'parameter_fitting':False, 'historical_investigations_rerun':False,
 'next_necessary_step':'Separately authorize reconsidering the centered far-field derivative restriction before any further feasibility study'
}
Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(out['exact_contradiction'])
