#!/usr/bin/env python3
"""Exact bounded-domain identities; no fitted values or numerical PDE solves."""
import sympy as S

r, m, R, ri, ro = S.symbols('r m R r_i r_o', positive=True)
F, g, gi, go = S.symbols('F g g_i g_o', real=True)
AI, AK = S.symbols('A_I A_K', real=True)
print('SymPy', S.__version__)


def reduce_bessel(expr):
    # Differentiate first, then apply the order-two recurrences exactly.
    expr = S.expand(expr.doit())
    for arg in (m*r, m*R, m*ri, m*ro):
        expr = expr.subs(S.besseli(2, arg), S.besseli(0, arg)-2*S.besseli(1, arg)/arg)
        expr = expr.subs(S.besselk(2, arg), S.besselk(0, arg)+2*S.besselk(1, arg)/arg)
    return S.simplify(expr)


def zero(label, expr):
    vals = list(expr) if isinstance(expr, S.MatrixBase) else [expr]
    vals = [reduce_bessel(e) for e in vals]
    assert all(e == 0 for e in vals), (label, vals)
    print(label, '=', vals)


def L(expr):
    return S.diff(expr, r, 2)+S.diff(expr, r)/r-m**2*expr


I, K = S.besseli(0, m*r), S.besselk(0, m*r)
p = -F*(r**2/m**2+4/m**4)
alpha, beta = S.symbols('alpha beta')
coeff = S.Poly(S.expand(L(alpha*r**2+beta)-F*r**2), r)
sol = S.solve(coeff.all_coeffs(), (alpha, beta), dict=True)
assert sol == [{alpha: -F/m**2, beta: -4*F/m**4}], sol
print('particular coefficients =', sol)
zero('particular residual', L(p)-F*r**2)
zero('I0 homogeneous residual', L(I))
zero('K0 homogeneous residual', L(K))
zero('I0 radial derivative', S.diff(I, r)-m*S.besseli(1, m*r))
zero('K0 radial derivative', S.diff(K, r)+m*S.besselk(1, m*r))
u = p+AI*I+AK*K
zero('general radial amplitude residual', L(u)-F*r**2)
# Disk regular mode is analytic in r^2; excludes logarithmic K0.
zero('I0 axis value', I.subs(r, 0)-1)
zero('I0 axis derivative', S.diff(I, r).subs(r, 0))
assert S.limit(K, r, 0, dir='+') == S.oo
print('K0 axis singularity = +oo')
uD = p+(g-p.subs(r, R))*I/I.subs(r, R)
uN = p+2*F*R*I/(m**3*S.besseli(1, m*R))
zero('disk Dirichlet trace', uD.subs(r, R)-g)
zero('disk Neumann derivative', S.diff(uN, r).subs(r, R))
zero('disk Dirichlet radial residual', L(uD)-F*r**2)
zero('disk Neumann radial residual', L(uN)-F*r**2)
t = m*R
zero('disk zero-trace derivative obstruction expression',
     S.diff(uD.subs(g, 0), r).subs(r, R)
     -F/m**3*((t**2+4)*S.besseli(1, t)/S.besseli(0, t)-2*t))
Dmat = S.Matrix([[I.subs(r, ri), K.subs(r, ri)],
                 [I.subs(r, ro), K.subs(r, ro)]])
dvec = S.Matrix([gi-p.subs(r, ri), go-p.subs(r, ro)])
cD = Dmat.inv()*dvec
uAD = u.subs({AI:cD[0], AK:cD[1]})
zero('annular Dirichlet coefficient system', Dmat*cD-dvec)
zero('annular inner Dirichlet trace', uAD.subs(r, ri)-gi)
zero('annular outer Dirichlet trace', uAD.subs(r, ro)-go)
Nmat = S.Matrix([[m*S.besseli(1, m*ri), -m*S.besselk(1, m*ri)],
                 [m*S.besseli(1, m*ro), -m*S.besselk(1, m*ro)]])
nvec = S.Matrix([2*F*ri/m**2, 2*F*ro/m**2])
cN = Nmat.inv()*nvec
uAN = u.subs({AI:cN[0], AK:cN[1]})
zero('annular Neumann coefficient system', Nmat*cN-nvec)
zero('annular inner outward Neumann', -S.diff(uAN, r).subs(r, ri))
zero('annular outer outward Neumann', S.diff(uAN, r).subs(r, ro))
# Avoid expanding huge coefficients: operator is linear and coefficients constant.
for label, cc in [('Dirichlet', cD), ('Neumann', cN)]:
    zero('annular '+label+' radial residual',
         L(p)-F*r**2+cc[0]*reduce_bessel(L(I))+cc[1]*reduce_bessel(L(K)))
zero('Green matching identity', S.diff(r*(I*S.diff(u, r)-u*S.diff(I, r)), r)-F*I*r**3)
primitive = r**3*S.besseli(1, m*r)/m-2*r**2*S.besseli(0, m*r)/m**2+4*r*S.besseli(1, m*r)/m**3
zero('positive-weight forcing integral primitive', S.diff(primitive, r)-I*r**3)
# Integrating (r u')'=m^2 r u+F r^3 gives the Neumann means.
zero('integrated amplitude identity', S.diff(r*S.diff(u, r), r)-m**2*r*u-F*r**3)
zero('disk Neumann mean formula', -F*(R**4/4)/(m**2*(R**2/2))+F*R**2/(2*m**2))
zero('annular Neumann mean formula', -F*(ro**4-ri**4)/(2*m**2*(ro**2-ri**2))+F*(ri**2+ro**2)/(2*m**2))

# Independent direct divergence of each linearized phase-current density.
# Use arbitrary transverse w; this includes every radial response above.
tt, x, y, z = X = S.symbols('t x y z', real=True)
rho, lam = S.symbols('rho_0 lambda_V', positive=True)
a, b, c, kt, kp, G, B, eps = S.symbols('a b c k_Theta k_Psi G B epsilon', real=True)
A = S.Matrix([[a, c], [c, b]])
M = S.eye(2)+rho**2*A
k = S.Matrix([kt, kp])
eta = S.diag(-1, 1, 1, 1)
s = x*x+y*y
# Inherit ordinary h from the verified baseline; no Einstein recalculation.
h = S.zeros(4)
h[0, 0] = 2*S.pi*G*(2*eps-B)*s
h[3, 3] = -2*S.pi*G*(B+2*eps)*s
for i in (1, 2):
    for j in (1, 2):
        h[i, j] = 2*S.pi*G*(B-eps)*s*int(i == j)-4*S.pi*G*eps*X[i]*X[j]
tr = 4*S.pi*G*(B-4*eps)*s
hu = eta*h*eta
Dz = sum(S.diff(hu[i, 3], X[i]) for i in range(4))-S.diff(tr, z)/2
zero('inherited metric D^z on response', Dz)
w = S.Function('w')(x, y)
f = S.Matrix([S.Function('f_Theta')(*X), S.Function('f_Psi')(*X)])
box = lambda e: -S.diff(e, tt, 2)+sum(S.diff(e, v, 2) for v in (x, y, z))
tau = S.symbols('tau')
Mt = S.eye(2)+(rho+tau*w)**2*A
inv = eta-tau*hu
root = 1+tau*tr/2
P = [k*int(i == 3)+tau*f.diff(X[i]) for i in range(4)]
phase = S.Matrix([S.diff(sum(S.diff(root*sum(inv[i,j]*(Mt*P[j])[n] for j in range(4)), X[i]) for i in range(4)), tau).subs(tau, 0) for n in range(2)])
expected = M*f.applyfunc(box)+2*rho*A*k*S.diff(w, z)-M*k*Dz
for n, label in enumerate(('Theta', 'Psi')):
    zero(label+' general phase perturbation equation', phase[n]-expected[n])
    zero(label+' fixed-phase current residual', phase[n].subs({f[0]:0, f[1]:0}).doit())
C = (k.T*A*k)[0]
forcing = 2*S.pi*G*rho*C*(B+2*eps)
zero('radial forcing agrees with inherited h_zz', -rho*C*h[3,3]-forcing*s)
# The radial solution has m^2=2 lambda rho_0^2, F=forcing.
zero('full first-order fixed-phase radial amplitude residual',
     (L(u)-F*r**2).subs({F:forcing, m:S.sqrt(2*lam*rho**2)}))
print('PASS: exact bounded response, boundary and both phase-current identities.')
print('Analytical positivity/coercivity and matching proofs are in the companion documents.')
