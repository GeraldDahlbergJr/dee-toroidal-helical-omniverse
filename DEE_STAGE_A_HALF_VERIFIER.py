#!/usr/bin/env python3
"""Exact local first-order checks; no numerical solving. Run with Python + SymPy."""
import sympy as S

print('SymPy', S.__version__)
t, x, y, z = X = S.symbols('t x y z', real=True)
G, eps, B = S.symbols('G epsilon B', real=True)
r, lam = S.symbols('rho_0 lambda_V', positive=True)
a, b, c, p, q = S.symbols('a b c k_Theta k_Psi', real=True)
eta = S.diag(-1, 1, 1, 1)
sg = [-1, 1, 1, 1]
A = S.Matrix([[a, c], [c, b]])
k = S.Matrix([p, q])
M = S.eye(2) + r**2*A
C = (k.T*A*k)[0]
B0 = (k.T*M*k)[0]
v2 = r**2 + C/lam  # exact equilibrium substitution, no square roots
V0 = lam*(r**2-v2)**2/4
eps0 = V0+B0/2

def zero(label, expr):
    vals = list(expr) if isinstance(expr, S.MatrixBase) else [expr]
    vals = [S.simplify(e) for e in vals]
    if any(e != 0 for e in vals):
        raise AssertionError((label, vals))
    print(label, '=', vals if len(vals) > 1 else vals[0])

def box(f):
    return sum(sg[i]*S.diff(f, X[i], 2) for i in range(4))

s = x*x+y*y
hb = S.zeros(4)
hb[0, 0] = -4*S.pi*G*eps*s
hb[3, 3] = -4*S.pi*G*(B-eps)*s
for i in (1, 2):
    for j in (1, 2):
        hb[i, j] = 6*S.pi*G*eps*s*int(i == j)-4*S.pi*G*eps*X[i]*X[j]
T0 = S.diag(eps, -eps, -eps, B-eps)
for i in range(4):
    for j in range(i, 4):
        zero('Einstein '+str((i, j)), box(hb[i, j])+16*S.pi*G*T0[i, j])
for j in range(4):
    zero('gauge '+str(j), sum(sg[i]*S.diff(hb[i, j], X[i]) for i in range(4)))
trb = S.trace(eta*hb)
h = S.simplify(hb-eta*trb/2)
zero('trace reversal involution', h-eta*S.trace(eta*h)/2-hb)
print('bar trace =', S.factor(trb))
print('h =', h.applyfunc(S.factor))
# Independent linearized Ricci/Einstein calculation from ordinary h.
trh = S.trace(eta*h)
Ric = S.Matrix(4, 4, lambda i,j: S.Rational(1,2)*(
    sum(sg[d]*(S.diff(h[d,j], X[d], X[i])+S.diff(h[d,i], X[d], X[j])) for d in range(4))
    - box(h[i,j])-S.diff(trh, X[i], X[j])))
zero('direct linearized Einstein tensor', Ric-eta*S.trace(eta*Ric)/2-8*S.pi*G*T0)
zero('background amplitude', -lam*r*(r*r-v2)-r*C)
zero('background stress from action', S.diag(0,0,0,B0)-eta*eps0-T0.subs({B:B0,eps:eps0}))

# Arbitrary spacetime-dependent metric and both matter perturbations.
H = S.zeros(4)
for i in range(4):
    for j in range(i,4):
        H[i,j] = H[j,i] = S.Function('h'+str(i)+str(j))(*X)
u = S.Function('delta_rho')(*X)
f = S.Matrix([S.Function('delta_Theta')(*X), S.Function('delta_Psi')(*X)])
tr = S.trace(eta*H)
Hu = eta*H*eta
D = [sum(S.diff(Hu[i,j], X[i]) for i in range(4))-sg[j]*S.diff(tr,X[j])/2 for j in range(4)]
# Coefficient extraction is exact differentiation at tau=0, not finite differences.
tau = S.symbols('tau')
R = r+tau*u
Mt = S.eye(2)+R**2*A
inv = eta-tau*Hu
root = 1+tau*tr/2
P = [k*int(i==3)+tau*f.diff(X[i]) for i in range(4)]
coeff = lambda e: S.simplify(S.diff(e,tau).subs(tau,0))
phase = S.Matrix([coeff(sum(S.diff(root*sum(inv[i,j]*(Mt*P[j])[n] for j in range(4)),X[i]) for i in range(4))) for n in range(2)])
phase_expected = M*f.applyfunc(box)+2*r*A*k*S.diff(u,z)-M*k*D[3]
zero('both general phase equations from current density', phase-phase_expected)
kinA = sum(inv[i,j]*(P[i].T*A*P[j])[0] for i in range(4) for j in range(4))
Er = coeff(sum(S.diff(root*sum(inv[i,j]*S.diff(R,X[j]) for j in range(4)),X[i]) for i in range(4))/root-lam*R*(R**2-v2)-R*kinA)
Er_expected = box(u)-2*lam*r*r*u+r*C*H[3,3]-2*r*(k.T*A*f.diff(z))[0]
zero('general amplitude equation from action', Er-Er_expected)
# Stress variation independently extracted from the action stress formula.
kinQ = sum(inv[i,j]*(P[i].T*Mt*P[j])[0] for i in range(4) for j in range(4))
Lm = -sum(inv[i,j]*S.diff(R,X[i])*S.diff(R,X[j]) for i in range(4) for j in range(4))/2-lam*(R**2-v2)**2/4-kinQ/2
DT = S.Matrix(4,4,lambda i,j: coeff(S.diff(R,X[i])*S.diff(R,X[j])+(P[i].T*Mt*P[j])[0]+(eta[i,j]+tau*H[i,j])*Lm))
d = (k.T*M*f.diff(z))[0]
DT_expected = S.Matrix(4,4,lambda i,j: 2*r*C*u*int(i==3)*int(j==3)+(k.T*M*f.diff(X[j]))[0]*int(i==3)+(k.T*M*f.diff(X[i]))[0]*int(j==3)-eps0*H[i,j]-eta[i,j]*(d-B0*H[3,3]/2))
zero('all stress variations from action', DT-DT_expected)
Tbase = T0.subs({B:B0,eps:eps0})
Tmix = eta*Tbase
DTmix = eta*DT-Hu*Tbase
Gamma = lambda l,i,j: sg[l]*(S.diff(H[l,j],X[i])+S.diff(H[l,i],X[j])-S.diff(H[i,j],X[l]))/2
cons = []
for n in range(4):
    div = sum(S.diff(DTmix[i,n],X[i]) for i in range(4))
    div += sum(Gamma(i,i,l)*Tmix[l,n]-Gamma(l,i,n)*Tmix[i,l] for i in range(4) for l in range(4))
    cons.append(S.simplify(div))
    zero('conservation identity '+str(n), div-(k.T*phase)[0]*int(n==3))
# Frozen matter: both phase equations and all conservation components still vanish.
subh = {H[i,j]:h[i,j] for i in range(4) for j in range(i,4)}
freeze = {u:0, f[0]:0, f[1]:0}
zero('frozen phase residuals', phase.subs(freeze).subs(subh).doit())
zero('frozen conservation residuals', S.Matrix(cons).subs(freeze).subs(subh).doit())
frozen = S.factor(Er.subs(freeze).subs(subh).doit())
zero('reported frozen amplitude residual', frozen+2*S.pi*G*r*C*(B+2*eps)*s)
print('frozen E_rho =', frozen)
# Verify the transverse reduction with an arbitrary response, not a PDE solve.
w = S.Function('transverse_response')(x,y)
transverse = {u:w, f[0]:0, f[1]:0}
zero('fixed-phase transverse amplitude reduction',
     Er.subs(transverse).subs(subh).doit()
     -(S.diff(w,x,2)+S.diff(w,y,2)-2*lam*r*r*w+r*C*h[3,3]))
zero('fixed-phase transverse phase equations',
     phase.subs(transverse).subs(subh).doit())
# RHS and residual conventions explicitly distinguished.
print('fixed-phase transverse RHS =', S.factor(-r*C*h[3,3]))
print('PASS: all exact symbolic checks; no global existence or localization claim.')
