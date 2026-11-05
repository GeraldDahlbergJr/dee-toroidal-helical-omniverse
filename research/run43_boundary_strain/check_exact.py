"""One SBP boundary strain candidate; exact bounded checks, no field solves."""
import json
import sympy as s
R = s.Rational
n = 7
D = s.zeros(n)
D[0,0], D[0,1] = -1, 1
D[-1,-2], D[-1,-1] = -1, 1
for i in range(1,n-1):
    D[i,i-1], D[i,i+1] = -R(1,2), R(1,2)
H = s.diag(R(1,2),1,1,1,1,1,R(1,2))
B = s.diag(-1,0,0,0,0,0,1)
assert H*D+D.T*H == B
I = s.eye(n)[:,1:-1]
K = (D*I).T*H*(D*I)
assert I.T*H*D*D*I == -K
assert K == K.T
minors = [K[:j,:j].det() for j in range(1,n-1)]
assert all(m > 0 for m in minors)
# Local actual-toroidal radial flux, no angular grid assembled.
h = s.symbols('h', positive=True)
a,b,z,r = s.symbols('a b z r', real=True)
u = (r-a)*(b-r)
J = r*(1+r*z)
up = s.diff(u,r)
F = J*up
# Endpoint forward/backward derivative errors on this quadratic: -h,+h.
def at(expr,x): return expr.subs(r,x)
inner = (at(F,a+2*h)-at(F,a)+h*at(J,a))/(2*h*at(J,a+h))
inner2 = (at(F,a+3*h)-at(F,a+h))/(2*h*at(J,a+2*h))
outer = (at(F,b)-at(F,b-2*h)+h*at(J,b))/(2*h*at(J,b-h))
outer2 = (at(F,b-h)-at(F,b-3*h))/(2*h*at(J,b-2*h))
exact = s.diff(F,r)/J
errors = [s.factor(v-at(exact,x)) for v,x in zip(
    [inner,inner2,outer,outer2],[a+h,a+2*h,b-h,b-2*h])]
limits = [s.simplify(s.limit(e,h,0)) for e in errors]
assert limits == [R(1,2),0,R(1,2),0]
second = [s.simplify(s.limit(errors[j]/h**2,h,0)) for j in (1,3)]
assert second == [-2*z/(a*(1+a*z)),-2*z/(b*(1+b*z))]
actual = [s.simplify(v*R(4,3)) for v in limits]
assert actual == [R(2,3),0,R(2,3),0]
result = dict(base_commit='b075780b0c30296f0e7fbd6fb2230d503e31b40e',
    candidate='minimal diagonal-norm SBP endpoint strain', arithmetic='exact symbolic/rational',
    sympy=s.__version__, radial_nodes=n, largest_matrix=[n,n],
    scalar_SBP_identity=True, Dirichlet_Gram_identity=True,
    normal_Gram_leading_minors=list(map(str,minors)),
    toroidal_weighted_adjoint='PROVEN algebraically in derivation.md',
    radial_errors=list(map(str,errors)),
    radial_error_limits_inner1_inner2_outer1_outer2=list(map(str,limits)),
    actual_x_momentum_error_limits=list(map(str,actual)),
    second_row_radial_error_over_h_squared=list(map(str,second)),
    boundary_closure='FAILED: nonzero O(1) first-interior-row error',
    coercivity='PARTIAL: normal reduction positive; full toroidal Korn bound unresolved',
    physical_solves=0, full_grid_sweeps=0,
    stop='Consistency failure: no further candidates or physical tests')
if __name__ == '__main__':
    print(json.dumps(result,indent=2))
