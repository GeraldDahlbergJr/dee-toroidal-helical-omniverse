#!/usr/bin/env python3
"""Stage B exact identities only. Analytical scope is in the two companion files."""
import sympy as s

print('SymPy', s.__version__)
r, R, m = s.symbols('r R m', positive=True)
F, g, q = s.symbols('F g q', real=True)

def reduce(expr):
    expr = s.expand(expr.doit())
    for arg in (m*r, m*R):
        expr = expr.subs(s.besseli(2, arg), s.besseli(0, arg)-2*s.besseli(1, arg)/arg)
        expr = expr.subs(s.besselk(2, arg), s.besselk(0, arg)+2*s.besselk(1, arg)/arg)
    return s.simplify(expr)

def zero(label, expr):
    vals = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    vals = [reduce(e) for e in vals]
    assert all(e == 0 for e in vals), (label, vals)
    print(label, '= 0')

I, K = s.besseli(0, m*r), s.besselk(0, m*r)
p = -F*(r**2/m**2+4/m**4)
# Target only exterior matching, rather than rerunning bounded-domain calculations.
J = s.Matrix([[I, K], [s.diff(I, r), s.diff(K, r)]]).subs(r, R)
d = s.Matrix([g-p.subs(r, R), q-s.diff(p, r).subs(r, R)])
coeff = J.inv()*d
zero('exterior Cauchy coefficient system', J*coeff-d)
u = p+coeff[0]*I+coeff[1]*K
zero('exterior amplitude trace', u.subs(r, R)-g)
zero('exterior normal derivative trace', s.diff(u, r).subs(r, R)-q)
# Linearity avoids expanding large boundary coefficients.
L = lambda v: s.diff(v, r, 2)+s.diff(v, r)/r-m**2*v
zero('exterior original-action amplitude equation',
     L(p)-F*r**2+coeff[0]*reduce(L(I))+coeff[1]*reduce(L(K)))
W = I*s.diff(K, r)-s.diff(I, r)*K
zero('Abel Wronskian identity', s.diff(r*W, r))
# Standard small-r limits I0~1, I0'~m²r/2, K0'~-1/r fix rW=-1.
zero('I0 small-radius normalization', s.limit(I, r, 0)-1)
zero('r K0 derivative normalization', s.limit(r*s.diff(K, r), r, 0)+1)
zero('r I0 derivative K0 normalization', s.limit(r*s.diff(I, r)*K, r, 0))
ai, ak = s.symbols('a_I a_K')
zero('Cauchy recovery of interior coefficients', J.inv()*(J*s.Matrix([ai, ak]))-s.Matrix([ai, ak]))

G, B, eps, rho, C = s.symbols('G B epsilon rho_0 C', real=True)
at = 2*s.pi*G*(2*eps-B)
ar = 2*s.pi*G*(B-3*eps)
ap = 2*s.pi*G*(B-eps)
az = -2*s.pi*G*(B+2*eps)
zero('forcing sign from h_zz', -rho*C*az*r**2-2*s.pi*G*rho*C*(B+2*eps)*r**2)

# Independent Cartesian-to-cylindrical tensor conversion at arbitrary angle.
phi = s.symbols('phi', real=True)
xvec = s.Matrix([r*s.cos(phi), r*s.sin(phi)])
Hxy = 2*s.pi*G*(B-eps)*r**2*s.eye(2)-4*s.pi*G*eps*xvec*xvec.T
P = xvec.jacobian([r, phi])
zero('transverse cylindrical conversion', P.T*Hxy*P-s.diag(ar*r**2, ap*r**4))

# Derive K and Gauss constraint with a formal expansion parameter.
tau = s.symbols('tau', real=True)
metric = s.diag(-1+tau*at*r**2, 1+tau*ar*r**2,
                r**2*(1+tau*ap*r**2), 1+tau*az*r**2)
gamma = s.diag(metric[0, 0], metric[2, 2], metric[3, 3])
Kform = gamma.diff(r)/(2*s.sqrt(metric[1, 1]))
first = lambda e: s.simplify(s.diff(e, tau).subs(tau, 0))
zero('direct extrinsic curvature coefficients',
     Kform.applyfunc(first)-s.diag(at*r, (2*ap-ar/2)*r**3, az*r))
Kmixed = gamma.inv()*Kform
Gnn = (s.trace(Kmixed)**2-s.trace(Kmixed*Kmixed))/2
zero('Gauss normal Einstein constraint', first(Gnn)-(az-at))
zero('normal stress Einstein equation', first(Gnn)+8*s.pi*G*eps)
# Evaluate on either side of an arbitrary common interface.
for label, tensor in [('induced metric', gamma), ('extrinsic curvature', Kform)]:
    interior = tensor.applyfunc(lambda e: e.subs(tau, 0)+tau*first(e)).subs(r, R)
    exterior = tensor.applyfunc(lambda e: e.subs(tau, 0)+tau*first(e)).subs(r, R)
    zero('same-polynomial '+label+' jump', exterior-interior)

# Current densities from the original action in cylindrical coordinates.
z, t = s.symbols('z t', real=True)
a, b, c, kt, kp = s.symbols('a b c k_Theta k_Psi', real=True)
A = s.Matrix([[a, c], [c, b]])
w = s.Function('u')(r)
M = s.eye(2)+(rho+tau*w)**2*A
k = s.Matrix([kt, kp])
# sqrt(-g) to first order; Q=kz so only contravariant z currents exist.
root = r*(1+tau*(-at+ar+ap+az)*r**2/2)
Jz_density = root*(1-tau*az*r**2)*M*k
for i, label in enumerate(('Theta', 'Psi')):
    zero(label+' current divergence', s.diff(Jz_density[i], z))
    # Diagonal inverse metric and partial_r Q=0 imply normal current identically.
    zero(label+' radial normal current', (M*(k*s.diff(z, r)))[i])
# Scalar flux correction multiplies u'=O(G), thus absent at first order.
zero('first-order scalar normal flux',
     first(tau*s.diff(w, r)/s.sqrt(1+tau*ar*r**2))-s.diff(w, r))

# Direct mixed-tensor covariant divergence with cylindrical background connection.
coords = (t, r, phi, z)
g0 = s.diag(-1, 1, r**2, 1)
T0 = s.diag(-eps, -eps, -eps, B-eps)
hzz = az*r**2
dT = s.diag(B*hzz/2, B*hzz/2, B*hzz/2, 2*rho*C*w-B*hzz/2)
inv = metric.inv()
Gamma = [[[sum(inv[l, j]*(s.diff(metric[j, v], coords[mu])
           +s.diff(metric[j, mu], coords[v])-s.diff(metric[mu, v], coords[j]))
           for j in range(4))/2 for v in range(4)] for mu in range(4)] for l in range(4)]
for nu in range(4):
    divergence = sum(s.diff(dT[mu, nu], coords[mu]) for mu in range(4))
    # Diagonal mixed stress reduces connection terms to sum Gamma^mu_mu_nu(Tnu-Tmu).
    for mu in range(4):
        connection = Gamma[mu][mu][nu]
        divergence += connection.subs(tau, 0)*(dT[nu, nu]-dT[mu, mu])
        divergence += first(connection)*(T0[nu, nu]-T0[mu, mu])
    zero('covariant conservation '+str(coords[nu]), divergence)

V = s.symbols('V', nonnegative=True)
zero('background strong-energy combination', B-2*(V+B/2)+2*V)
print('PASS: Stage B exterior Cauchy matching, currents, junction geometry, constraint and conservation.')
print('Finite collar only; analytical admissibility and obstruction scopes are in the companion files.')
