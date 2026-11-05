"""Exact moment contradiction only; no closure implementation or physical solves."""
import json
import sympy as s

t = s.symbols('t')
R = s.Rational
# Local symbolic rows needed to evaluate composition at rows 1 through 4.
d = s.zeros(6, 7)
d[:4, :5] = s.Matrix([
    [-1, 2*t, 2-4*t, 2*t-1, 0],
    [-t, 0, 3*t-1, 1-2*t, 0],
    [2*t-1, 1-3*t, 0, t, 0],
    [R(1, 2)-t, 2*t-1, -t, 0, R(1, 2)]])
for i in (4, 5):
    d[i, i-1], d[i, i+1] = -R(1, 2), R(1, 2)
assert d*s.ones(7, 1) == s.zeros(6, 1)
assert d*s.Matrix(range(7)) == s.ones(6, 1)

def residual(k, i):
    first = d*s.Matrix([j**k for j in range(7)])
    value = sum(d[i, j]*first[j] for j in range(6))
    target = k*(k-1)*i**(k-2) if k >= 2 else 0
    return s.factor(value-target)

quadratic = [residual(2, i) for i in (1, 2)]
assert all(s.expand(a-b) == 0 for a, b in zip(
    quadratic, [2*(3*t-2)*(3*t-1), 4*(2*t-1)*(3*t-1)]))
gcd = s.monic(s.gcd(*quadratic), t)
assert gcd == t-R(1, 3)
forced = R(1, 3)
transition = [residual(2, i).subs(t, forced) for i in (3, 4)]
assert transition == [R(1, 18), -R(1, 6)]
cubic = [residual(3, i).subs(t, forced) for i in (1, 2)]
assert cubic == [R(14, 3), -R(4, 3)]
result = dict(
    base_commit='a6334b7c5757e56748de6b23cd55101e1f9df317',
    candidate='four-row trapezoidal-norm SBP strain block',
    arithmetic='exact rational/symbolic', sympy=s.__version__,
    quadratic_residuals_first_two=list(map(str, quadratic)),
    common_root=str(forced),
    quadratic_residuals_rows_3_4=list(map(str, transition)),
    cubic_residuals_first_two=list(map(str, cubic)),
    contradiction='first two rows force t=1/3; row 3 then requires 1/18=0',
    moment_conditions='FAIL', boundary_closure='FAILED',
    adjoint_identity='PASS: formal weighted identity proved in derivation',
    boundary_consistency='FAILED: nonzero O(1) error at row 3',
    coercivity='strict bound disproven by explicit Dirichlet strain kernel',
    physical_solves=0, grid_sweeps=0,
    stop='No further candidate or physical tests')
print(json.dumps(result, indent=2))
