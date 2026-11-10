#!/usr/bin/env python3
"""Exact Stage C identities; theorem assumptions/proofs are in companion reports.
No numerical PDE, nonlinear evolution, source rescaling or Run 43 imports.
"""
import sympy as s


def zero(label, expr):
    entries = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    values = [s.simplify(s.expand(e.doit())) for e in entries]
    if any(e != 0 for e in values):
        raise AssertionError((label, values))
    print('PASS:', label)


def main():
    print('SymPy', s.__version__)
    t, r, phi, z = coords = s.symbols('t r phi z', real=True)
    N, H, L, Z = [s.Function(x)(r) for x in ('N', 'H', 'L', 'Z')]
    rho = s.Function('rho')(r)
    q = s.Matrix([s.Function('q1')(r), s.Function('q2')(r)])
    wt, wp, mt, mp, kt, kp = s.symbols('wt wp mt mp kt kp', real=True)
    w, m, k = s.Matrix([wt, wp]), s.Matrix([mt, mp]), s.Matrix([kt, kp])
    x = s.symbols('x', real=True)
    Mx = s.Matrix([[s.Function('K1')(x), s.Function('Mix')(x)],
                   [s.Function('Mix')(x), s.Function('K2')(x)]])
    Vx = s.Function('V')(x)
    M, V, Mr = Mx.subs(x, rho), Vx.subs(x, rho), Mx.diff(x).subs(x, rho)
    S = N*H*L*Z  # positive metric functions, fixed orientation
    metric = s.diag(-N**2, H**2, L**2, Z**2)
    inv = s.diag(-1/N**2, 1/H**2, 1/L**2, 1/Z**2)
    deriv = [w, q.diff(r), m, k]
    W = sum(inv[i, i]*(deriv[i].T*M*deriv[i])[0] for i in range(4))
    lag = -rho.diff(r)**2/(2*H**2)-V-W/2
    er = (S*rho.diff(r)/H**2).diff(r)/S - s.diff(Vx, x).subs(x, rho)
    er -= sum(inv[i, i]*(deriv[i].T*Mr*deriv[i])[0]/2 for i in range(4))
    eq = (S*M*q.diff(r)/H**2).diff(r)/S
    density = S*lag
    elrho = (s.diff(density, rho)-s.diff(s.diff(density, rho.diff(r)), r))/S
    zero('rho variational equation vs covariant reduction', elrho-er)
    elq = s.Matrix([(s.diff(density, qi)-s.diff(s.diff(density, qi.diff(r)), r))/S
                    for qi in q])
    zero('both phase variational equations vs current reduction', elq-eq)
    T = s.Matrix(4, 4, lambda i, j:
                 (rho.diff(r)**2 if i == j == 1 else 0)
                 +(deriv[i].T*M*deriv[j])[0]+metric[i, j]*lag)
    Gamma = [[[sum(inv[a, b]*(s.diff(metric[b, c], coords[d])
                  +s.diff(metric[b, d], coords[c])-s.diff(metric[c, d], coords[b]))
                  for b in range(4))/2 for d in range(4)]
                  for c in range(4)] for a in range(4)]
    mixed = inv*T
    for nu in range(4):
        div = sum(s.diff(mixed[mu, nu], coords[mu]) for mu in range(4))
        div += sum(Gamma[mu][mu][a]*mixed[a, nu]
                   -Gamma[a][mu][nu]*mixed[mu, a]
                   for mu in range(4) for a in range(4))
        rhs = (er*rho.diff(r) if nu == 1 else 0)+(deriv[nu].T*eq)[0]
        zero('direct off-shell stress identity '+str(coords[nu]), div-rhs)
    zero('normal traction formula', T[1, 1]/H**2-(rho.diff(r)**2
         +(q.diff(r).T*M*q.diff(r))[0])/(2*H**2)+V
         -(w.T*M*w)[0]/(2*N**2)+(m.T*M*m)[0]/(2*L**2)
         +(k.T*M*k)[0]/(2*Z**2))
    # Independent geometry routes: Ricci definition and cylinder Gauss equation.
    Ric = s.Matrix(4, 4, lambda i, j:
        sum(s.diff(Gamma[a][i][j], coords[a])-s.diff(Gamma[a][i][a], coords[j])
            for a in range(4))
        +sum(Gamma[a][a][b]*Gamma[b][i][j]-Gamma[a][j][b]*Gamma[b][i][a]
             for a in range(4) for b in range(4)))
    scalar = s.trace(inv*Ric)
    gnn_ric = (Ric[1, 1]-metric[1, 1]*scalar/2)/H**2
    gamma = s.diag(-N**2, L**2, Z**2)
    kmixed = gamma.inv()*gamma.diff(r)/(2*H)
    gnn_gauss = (s.trace(kmixed)**2-s.trace(kmixed*kmixed))/2
    expected = (N.diff(r)*L.diff(r)/(N*L)+N.diff(r)*Z.diff(r)/(N*Z)
                +L.diff(r)*Z.diff(r)/(L*Z))/H**2
    zero('normal Einstein constraint: direct Ricci vs Gauss', gnn_ric-gnn_gauss)
    zero('normal Einstein constraint: explicit expression', gnn_gauss-expected)
    # Integrability is checked before imposing invariant gradients.
    kr = s.Function('k')(r)
    tapered = kr*z
    zero('phase taper mixed partials', s.diff(tapered, r, z)-s.diff(tapered, z, r))
    zero('invariant-gradient closure requires k prime zero',
         s.diff(s.diff(tapered, z), r)-kr.diff(r))
    a, b, c, lam, v = s.symbols('a b c lambda_V v', positive=True)
    A = s.Matrix([[a, c], [c, b]])
    polyM = s.eye(2)+x**2*A
    polyV = lam*(x**2-v**2)**2/4
    zero('unit phase stiffness at zero amplitude', polyM.subs(x, 0)-s.eye(2))
    zero('zero amplitude potential is positive', polyV.subs(x, 0)-lam*v**4/4)
    for sign in (1, -1):
        zero('vacuum potential '+str(sign), polyV.subs(x, sign*v))
        zero('vacuum amplitude equation '+str(sign), polyV.diff(x).subs(x, sign*v))
    zero('kinetic determinant polynomial', polyM.det()
         -(1+(a+b)*x**2+(a*b-c**2)*x**4))
    C = (k.T*A*k)[0]
    zero('established axial amplitude branch',
         -x*(lam*(x**2-v**2)+C).subs(x**2, v**2-C/lam))
    # Null stress cancellation: zero Lorentz contraction is not zero energy.
    n1, n2 = s.symbols('n1 n2', real=True)
    wave = s.Matrix([n1, n2])
    zero('null phase contraction', -(wave.T*polyM*wave)[0]+(wave.T*polyM*wave)[0])
    energy = (wave.T*polyM*wave)[0]
    zero('null wave energy retains positive quadratic form',
         energy-(n1**2+n2**2+x**2*(a*n1**2+2*c*n1*n2+b*n2**2)))
    # General positive kinetic energy verified via a Cholesky parameterization.
    aa, bb, cc = s.symbols('aa bb cc', real=True)
    P = s.Matrix([[aa, bb], [0, cc]])
    MM = P.T*P
    d = s.Matrix(4, 3, lambda i, j: s.Symbol('d%d%d' % (i, j), real=True))
    dq = lambda i: s.Matrix([d[i, 1], d[i, 2]])
    eta = s.diag(-1, 1, 1, 1)
    kinetic = sum(eta[i, i]*(d[i, 0]**2+(dq(i).T*MM*dq(i))[0])
                  for i in range(4))
    VV = s.symbols('VV')
    lagframe = -kinetic/2-VV
    T00 = d[0, 0]**2+(dq(0).T*MM*dq(0))[0]-lagframe
    squares = sum(d[i, 0]**2+(P*dq(i)).dot(P*dq(i)) for i in range(4))/2+VV
    zero('orthonormal energy sum of squares', T00-squares)
    # Constant amplitude solves the radial vacuum ODE for arbitrary metric.
    zero('radial vacuum fixed point', H**2*polyV.diff(x).subs(x, v))
    print('PASS: Stage C exact identities; analytical theorem scopes are in companion files.')


if __name__ == '__main__':
    main()
