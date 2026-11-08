# DEE matter independent check

Starting commit: 27ae206a7b25b2ba09fdc89cfe1bdbe7cc60dd6f.
Research date: 2026-11-07 America/Los_Angeles (2026-11-08 UTC).

## Verification methods and limits

The main derivation varies the covariant action by hand and reduces in cylindrical coordinates. This check instead uses (1) direct inverse-metric/field-gradient differentiation of the unreduced Cartesian Lagrangian, (2) constant Cartesian gradients for the exact solutions, and (3) explicit Christoffel contraction for all four off-shell cylindrical divergences. It does not call a solver, reuse a solver's derivative routine, or import repository scientific code. Exact symbolic simplification is not numerical residual certification.

Both derivation and checks were prepared by the same AI investigator. Independence here is independence of mathematical route and symbolic implementation, not independent personnel, peer review, or external scientific replication. Source transcription/common-model errors remain possible and warrant human review. The reproducible code is embedded below so a reviewer can run it separately. Environment used: Python 3 with SymPy 1.14.0. No PDE solve, numerical integration, parameter fitting, grid/refinement study or numerical certification was executed.

## A. Independent Euler–Lagrange and Hilbert checks in Cartesian coordinates

Treat the ten symmetric inverse-metric entries as independent variables, before setting g^mu_nu=diag(-1,1,1,1). The verifier builds L directly from Eq. (9). For a diagonal entry, -2 dL/dg^ii plus g_ii L is T_ii. For a single offdiagonal variable representing both g^ij and g^ji, -dL/dg^ij plus g_ij L gives T_ij. The distinct coefficient is necessary and is checked for all ten entries.

Differentiating L with respect to each phase gradient independently recovers both four-component currents with the manuscript's sign convention. Differentiating L with respect to rho recovers the amplitude force, while the canonical rho kinetic term supplies box rho upon integration by parts. Therefore the equations are checked against the original action rather than against the reduced formula alone.

For the axial exact solution, in Cartesian coordinates (t,x,y,z) the field gradients are constant: d rho=0, dTheta=(0,0,0,p), dPsi=(0,0,0,q). Every derivative of the conjugate gradient momenta vanishes. The entire amplitude Euler–Lagrange residual reduces to

$$-\rho[\lambda_V(\rho^2-v^2)+ap^2+2c_{\rm mix}pq+bq^2].$$

The verifier substitutes rho²=v²-(ap²+2c_mix pq+bq²)/lambda_V and obtains exactly zero. Cartesian phase momenta are constant, proving both phase equations directly without cylindrical first integrals. Independently differentiated T has exactly the diagonal form in the existence analysis and is constant, proving all four divergences vanish. Transforming x,y to r,phi yields T_rr=-epsilon, T_phiphi=-r²epsilon and zero mixed rphi: the r² factor is a coordinate transformation, not a physical singularity.

The null family is evaluated separately with gradients (p,0,0,p) and (q,0,0,q), rho=v. The amplitude residual vanishes and metric-differentiated stress equals N ell ell^T. Constant phase momenta and constant Cartesian stress independently prove its current and stress conservation.

## B. General off-shell cylindrical conservation check

The verifier treats rho(r), theta'(r), psi'(r) as arbitrary functions and each kinetic coefficient and V as arbitrary functions of rho, not only polynomials. It constructs stress, raises its first index with diag(-1,1,r^-2,1), and evaluates

$$\partial_\mu T^\mu{}_\nu+
\Gamma^\mu_{\mu\alpha}T^\alpha{}_\nu-
\Gamma^\alpha_{\mu\nu}T^\mu{}_\alpha.$$

For each nu, it subtracts E_rho partial_nu rho + E_Q^T partial_nu Q and symbolically obtains zero without imposing any field equation. In the radial component this independently confirms the sign and factors of
T_rr'+T_rr/r-T_phiphi/r³. In each tangential component it checks the nonzero radial phase-flux terms, not just the zero-current special case.

A second variational radial check can be made directly from the one-dimensional density
r L=-r rho'^2/2-r V+r w^TMw/2-r u^TMu/2-m^TMm/(2r)-r k^TMk/2.
Its Euler–Lagrange equations give r E_rho=0 and (r M u)'=0.
The unreduced variation is authoritative; this reduced-action check is used only for the retained scalar degrees of freedom.

## C. Non-computational checks and obstructions

Positivity does not rely on a sampled determinant: A positive definite implies every eigenvalue of I+rho²A exceeds one. Choosing k=s e with the explicit open bound in the existence analysis proves the allowed set is nonempty and rho_0 positive. The local winding theorem follows from a locally Lipschitz ODE right-hand side for r>0, rho>0 and det M>0; this is distinct from an exact global solution.

The constant-amplitude winding obstruction is independently visible in Cartesian gradients: |grad phi|²=1/(x²+y²); its amplitude-force coefficient is m^T A m. A constant V'(rho_0) and constant axial/time contractions cannot cancel that term on an open radial interval when m^T A m>0. The axis obstruction is also topological, independently of the energy-divergence argument. These checks do not assert a no-go theorem for variable-amplitude annular boundary-value problems.

## D. Reproduce the exact symbolic checks

Extract the following block to a temporary check.py and run python check.py in an environment with SymPy 1.14.0. It executes only the displayed algebra. No repository solver or archived certification is needed.

~~~python
import sympy as s
print('SymPy', s.__version__, '; exact symbolic algebra only')
R, lam, v, a, b, c, p, q = s.symbols('R lam v a b c p q', real=True)
M = s.Matrix([[1+a*R**2,c*R**2],[c*R**2,1+b*R**2]])
# Unreduced Cartesian matter Lagrangian; independent metric variation.
gvars = s.symbols('g00 g11 g22 g33 g01 g02 g03 g12 g13 g23')
g = s.Matrix([[gvars[0],gvars[4],gvars[5],gvars[6]],
 [gvars[4],gvars[1],gvars[7],gvars[8]],
 [gvars[5],gvars[7],gvars[2],gvars[9]],
 [gvars[6],gvars[8],gvars[9],gvars[3]]]) # inverse metric variables
flat = dict(zip(gvars,[-1,1,1,1,0,0,0,0,0,0]))
dr = s.Matrix(s.symbols('dR0:4')); dt = s.Matrix(s.symbols('dT0:4')); dp = s.Matrix(s.symbols('dP0:4'))
V = lam*(R**2-v**2)**2/4
L = -(dr.T*g*dr)[0]/2 - V - M[0,0]*(dt.T*g*dt)[0]/2 - M[1,1]*(dp.T*g*dp)[0]/2 - M[0,1]*(dt.T*g*dp)[0]
T = s.zeros(4)
pairs=[(0,0),(1,1),(2,2),(3,3),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
eta=s.diag(-1,1,1,1)
for gv,(i,j) in zip(gvars,pairs):
    # Offdiagonal variable changes both symmetric inverse metric entries.
    T[i,j]=T[j,i]=s.expand(((-2 if i==j else -1)*s.diff(L,gv)+eta[i,j]*L).subs(flat))
expected=dr*dr.T+M[0,0]*dt*dt.T+M[1,1]*dp*dp.T+M[0,1]*(dt*dp.T+dp*dt.T)+eta*L.subs(flat)
assert all(s.simplify(x)==0 for x in T-expected)
print('ten Hilbert stress components: verified by inverse-metric differentiation')
for i in range(4):
    assert s.simplify((-s.diff(L,dt[i])-(M[0,0]*g*dt+M[0,1]*g*dp)[i]).subs(flat))==0
    assert s.simplify((-s.diff(L,dp[i])-(M[1,1]*g*dp+M[0,1]*g*dt)[i]).subs(flat))==0
print('both four-component phase currents: verified by gradient differentiation')
C=a*p*p+2*c*p*q+b*q*q
branch={R**2:v**2-C/lam}
axial=dict(zip(dr,[0]*4)); axial.update(zip(dt,[0,0,0,p])); axial.update(zip(dp,[0,0,0,q]))
er=s.factor(s.diff(L,R).subs(flat).subs(axial))
assert s.simplify((er/R).subs(branch))==0
tax=T.subs(axial)
B=(s.Matrix([p,q]).T*M*s.Matrix([p,q]))[0]
eps=V+B/2
assert all(s.simplify(x)==0 for x in tax-s.diag(eps,-eps,-eps,B-eps))
print('axial exact family: amplitude equation and all Cartesian stress components verified')
null=dict(zip(dr,[0]*4)); null.update(zip(dt,[p,0,0,p])); null.update(zip(dp,[q,0,0,q]))
assert s.simplify(s.diff(L,R).subs(flat).subs(null).subs(R,v))==0
ell=s.Matrix([1,0,0,1]); N=(s.Matrix([p,q]).T*M.subs(R,v)*s.Matrix([p,q]))[0]
assert all(s.simplify(x)==0 for x in T.subs(null).subs(R,v)-N*ell*ell.T)
print('null axial exact family: amplitude equation and null-dust stress verified')
# Independent cylindrical connection check of all four stress divergences.
r=s.symbols('r',positive=True); rho=s.Function('rho')(r)
u=s.Matrix([s.Function('u')(r),s.Function('w')(r)])
om=s.Matrix(s.symbols('o0:2')); mm=s.Matrix(s.symbols('m0:2')); kk=s.Matrix(s.symbols('k0:2'))
H=s.Matrix([[s.Function('A')(rho),s.Function('C')(rho)],[s.Function('C')(rho),s.Function('B')(rho)]])
pot=s.Function('V')(rho)
bil=lambda x,y:(x.T*H*y)[0]
W=-bil(om,om)+bil(u,u)+bil(mm,mm)/r**2+bil(kk,kk)
lc=-s.diff(rho,r)**2/2-pot-W/2
vec=[om,u,mm,kk]; gc=s.diag(-1,1,r*r,1); inv=gc.inv()
tc=s.Matrix(4,4,lambda i,j:bil(vec[i],vec[j])+gc[i,j]*lc+(s.diff(rho,r)**2 if i==j==1 else 0))
# Explicit field-space contraction, independent of the reduced radial equation.
D=-om*om.T+u*u.T+mm*mm.T/r**2+kk*kk.T
er=s.diff(rho,r,2)+s.diff(rho,r)/r-s.diff(pot,rho)-sum(s.diff(H[i,j],rho)*D[i,j] for i in range(2) for j in range(2))/2
ep=(r*H*u).diff(r)/r
Gamma={(1,2,2):-r,(2,1,2):1/r,(2,2,1):1/r}
mixed=inv*tc
for nu in range(4):
    div=s.diff(mixed[1,nu],r)
    for mu in range(4):
        for al in range(4):
            div+=Gamma.get((mu,mu,al),0)*mixed[al,nu]-Gamma.get((al,mu,nu),0)*mixed[mu,al]
    rhs=(s.diff(rho,r)*er if nu==1 else 0)+(vec[nu].T*ep)[0]
    assert s.simplify(s.expand(div-rhs))==0
print('all four off-shell cylindrical Noether identities: verified using Christoffel contraction')
print('No PDE solves, numerical integration, or numerical certification performed.')
~~~


Observed output (all assertions completed):

~~~text
SymPy 1.14.0 ; exact symbolic algebra only
ten Hilbert stress components: verified by inverse-metric differentiation
both four-component phase currents: verified by gradient differentiation
axial exact family: amplitude equation and all Cartesian stress components verified
null axial exact family: amplitude equation and null-dust stress verified
all four off-shell cylindrical Noether identities: verified using Christoffel contraction
No PDE solves, numerical integration, or numerical certification performed.
~~~

## E. Preservation and final scope

The Git commit for this investigation adds only DEE_MATTER_FIELD_DERIVATION.md, DEE_MATTER_EXISTENCE_ANALYSIS.md and DEE_MATTER_INDEPENDENT_CHECK.md at repository root. Its parent is the exact requested starting commit. The original tree is reused and no pre-existing path is edited. This preserves manuscripts, the formal Eq. (34) correction, original audit findings, evidence seals, and every frozen Run 43 blob. The script is embedded in this document to keep the evidence package restricted to the three requested files.

Verdict: **EXISTENCE ESTABLISHED**, for the explicitly delimited nontrivial, extended, finite-local-energy matter configuration. Neither this symbolic check nor the existence theorem establishes a solution of Einstein's equations. No gravitational field is constructed or solved. Stop at this evidence commit; further gravitational work requires separate approval.
