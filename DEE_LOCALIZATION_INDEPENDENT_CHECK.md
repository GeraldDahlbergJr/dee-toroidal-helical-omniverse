# DEE localization independent mathematical checks

Starting commit: 14b8a1384d15c5e9fc296c2daf0f8d750b4dbee1.
Date: 2026-11-07 America/Los_Angeles (2026-11-08 UTC).

## 1. Meaning of independence

The proofs and checks were prepared by the same AI investigator. They use independent mathematical routes, not independent personnel or external peer review. Symbolic verification checks exact algebra, not analytic hypotheses or global existence. Each analytic hypothesis is stated in DEE_LOCALIZATION_ANALYSIS.md and DEE_BOUNDARY_FEASIBILITY.md. No numerical parameter selection, PDE solve, numerical integration, numerical certification, or gravitational construction is involved.

## 2. Claim-by-claim cross-check

| Claim | Primary method | Separate check and limits |
|---|---|---|
| w,k produce area-divergent energy | Uniform tail eigenvalue lower bound for M | Diagonalize the finite positive limiting M; each nonzero gradient has a positive limiting quadratic energy density. Its area integral diverges. Null Lorentzian contractions cannot cancel observer energy. |
| Winding has logarithmic tail/axis cost | Pointwise e>=delta |m|²/(2r²) | On each circle, Cauchy–Schwarz with integral d_phi Q=2pi m bounds integrated angular derivative energy below by pi delta |m|²/r. Integration in r diverges at the applicable endpoint. This applies also to nonuniform angular representatives with the same degree and uniform kinetic bound. |
| C must vanish | First integral and inverse-metric energy bound | Integrating each divergence equation over a disk gives constant radial flux 2pi C; a bounded source-free current through shrinking circles has zero flux. On an exterior tail the upper bound M<=L I gives the distinct logarithmic C-energy obstruction. |
| Only vacuum in the global ordinary-vacuum sector | Radial stress identity, with finite-energy subsequence | Cartesian transverse trace and cutoff dilation/Pohozaev give the same nonnegative potential integral; the first integrals then remove radial phase gradients and radial balance removes rho'. Two-dimensional gradient energy alone is scale invariant, so no unsupported universal sigma-model no-go is inferred. |
| No positive polynomial amplitude defect on a homogeneous transport background | Maximum/minimum principle on complete plane | The integrated radial ODE gives a monotonic departure argument from the central value; a positive-weight multiplier identity gives another check when its boundary term vanishes. These are detailed below. |
| No nontrivial smooth finite-radius vacuum patch | Source-free matching plus radial IVP uniqueness | Integrate field equations across a pillbox to determine derivative/current jumps; compare inward solution to the exact equilibrium solution with the same Cauchy data. Non-Lipschitz potentials and explicit interface sources are excluded. |
| Nonzero winding cannot fill a regular disk | Degree obstruction for circle-valued phases | Independent uniform-positive-M gradient energy lower bound diverges near the axis. An annulus avoids the topological contraction but supplies neither a regular core nor an integrable ordinary-vacuum tail. |
| Surface radial pressure jump | Distributional cylindrical connection calculation | Balance forces on a cylindrical surface element: interior positive pressure needs circumferential tension p_s=-R p_inside for vacuum outside, agreeing with [p_r]=p_s/R in the adopted orientation. |
| Annular local feasibility is preserved | Existing regular-radius ODE local-existence result | Explicit same-family interior/annulus/exterior continuation verifies all traces and fluxes but has infinite E_perp. This is not an existence proof for arbitrary two-boundary data. |

The area lower bound is
E>=pi delta (|w|²+|k|²)(R²-r_*²)/2; logarithmic coefficients are
pi delta |m|² and pi |C|²/L respectively. These are direct integrals with fixed lower bounds, not asymptotic fit coefficients.

### Cartesian Pohozaev check and boundary control

For static transverse scalars with target metric diag(1,M), write the full gradient Gram tensor H_ij=partial_i rho partial_j rho+partial_i Q^TM partial_j Q and L=-tr(H)/2-U. Metric variation gives T_ij=H_ij+delta_ij L; therefore
T_xx+T_yy=-2U exactly. Conservation and
partial_i(x^j T_ij)=T_xx+T_yy give, on a disk,

$$2pi R^2 p_r(R)=-2\int_{\text{disk}} U\,d^2x.$$

This is the radial identity obtained in another basis via a Cartesian divergence theorem. Finite nonnegative energy yields a sequence of vanishing boundary flux by |p_r|<=e and liminf R²e(R)=0. Thus no unproved universal falloff rate is needed. Smoothness at the origin removes the inner boundary contribution. Once U=0, radial balance gives r²p_r=constant=0 and the already reduced scalar sector has rho'=0. Annuli retain both boundary fluxes and cannot be treated as disks.

### Homogeneous-background no-defect theorem: independent ODE check

With positive rho_b, the polynomial equation is
(r rho')'=r lambda_V rho(rho²-rho_b²).
If rho_c>rho_b, its RHS starts positive and rho' becomes positive; as long as rho>rho_b this derivative sign cannot reverse because r rho' increases from zero. Hence rho never approaches the lower value rho_b. If 0<rho_c<rho_b, rho' becomes and remains negative while 0<rho<rho_b. It cannot approach the larger rho_b without crossing a turning point whose derivative sign the equation forbids; a crossing of zero would violate strict positive amplitude. If rho_c=rho_b, uniqueness of the regular-axis integral equation gives the constant solution: for a locally Lipschitz force, the difference of two solutions on [0,h] is bounded by (L h²/4) times its sup norm, forcing equality for small h and then regular-radius uniqueness extends it. These cases prove the same result as the maximum-principle argument without assuming derivative decay at infinity.

A further exact multiplier identity is

$$[r\rho'(\rho-\rho_b)]'
=r\rho'^2+r\lambda_V\rho(\rho-\rho_b)^2(\rho+\rho_b)
+r(\rho-\rho_b)E_{\rm bg},$$

where E_bg=rho''+rho'/r-lambda_V rho(rho²-rho_b²).
On shell with a vanishing boundary term its integral has two nonnegative terms for rho,rho_b>0, forcing rho=rho_b. This identity is not used to infer a boundary limit without proof; the ODE/maximum-principle proof covers the stated bounded-asymptotic case directly.

### Shell conditions and source-free uniqueness

The normal traction delta coefficient is [p_r]+Gamma^r_phi_phi S^{phi phi}
=[p_r]-p_s/R. Surface energy density and axial pressure do not enter this static normal curvature force, but still require a physically specified layer if a shell is proposed. Equal normal traction is necessary and not sufficient: continuous phase traces require matching w,k,m. Integrating the original scalar equations across an interface verifies flux matching independently of traction. Ordinary C², invertible-metric radial IVP uniqueness then forbids compact vacuum patches.

## 3. Reproducible symbolic checks

The following stand-alone code builds the radial action directly, derives its Euler–Lagrange amplitude residual, and compares it to the stress force balance. A separate Cartesian Gram-tensor calculation checks the transverse trace without importing the radial stress formulas. Further checks cover the effective potential, signed relative energy, multiplier identity and shell sign.

Run with Python 3 and SymPy 1.14.0. This is algebra only; it imports no repository solver. The script is embedded to keep the committed package to the four requested documents.

~~~python
import sympy as s
print("SymPy", s.__version__, "- exact localization identities")
r=s.symbols("r", positive=True)
rho=s.Function("rho")(r)
lam,v,a,b,c,rb=s.symbols("lam v a b c rb", positive=True)
A=s.Matrix([[a,c],[c,b]])
M=s.eye(2)+rho**2*A
w=s.Matrix(s.symbols("w0:2")); m=s.Matrix(s.symbols("m0:2"))
k=s.Matrix(s.symbols("k0:2"))
u=s.Matrix([s.Function("u0")(r),s.Function("u1")(r)])
bil=lambda x,y:(x.T*M*y)[0]
V=lam*(rho**2-v**2)**2/4
L=-s.diff(rho,r)**2/2-V+(bil(w,w)-bil(u,u)-bil(m,m)/r**2-bil(k,k))/2
pr=s.diff(rho,r)**2+bil(u,u)+L
pp=bil(m,m)/r**2+L
# Derive E_rho anew from radial action, not from stress balance.
radial_L=r*L
E=(s.diff(radial_L,rho)-s.diff(s.diff(radial_L,s.diff(rho,r)),r))/r
EQ=(r*M*u).diff(r)/r
assert s.simplify(s.expand(s.diff(pr,r)+(pr-pp)/r-s.diff(rho,r)*E-(u.T*EQ)[0]))==0
assert s.simplify(pr+pp-(bil(w,w)-bil(k,k)-2*V))==0
assert s.simplify(s.expand(s.diff(r*r*pr,r)-r*(bil(w,w)-bil(k,k)-2*V)-r*r*(s.diff(rho,r)*E+(u.T*EQ)[0])))==0
print("general cylindrical virial and off-shell force balance: verified")
# Independent Cartesian trace of static stress in two transverse dimensions.
hxx,hxy,hyy,U=s.symbols("hxx hxy hyy U")
# h_ij is the full scalar target-metric gradient Gram tensor.
Ld=-(hxx+hyy)/2-U
Txx=hxx+Ld; Tyy=hyy+Ld
assert s.expand(Txx+Tyy)==-2*U
print("independent Cartesian transverse stress trace: verified")
# Two-dimensional dilation of gradient and potential energies.
scale,Eg,Ep=s.symbols("scale Eg Ep", positive=True)
assert s.diff(Eg+scale**2*Ep,scale).subs(scale,1)==2*Ep
print("Derrick derivative: verified")
x=s.symbols("x", real=True)
d=lam*(v*v-rb*rb)
Ueff=lam*(x*x-v*v)**2/4+d*x*x/2
assert s.factor(Ueff-Ueff.subs(x,rb)-lam*(x*x-rb*rb)**2/4)==0
assert s.factor(s.diff(Ueff,x)-lam*x*(x*x-rb*rb))==0
S,W=s.symbols("S W")
e_rel=(lam*(x*x-v*v)**2/4-lam*(rb*rb-v*v)**2/4)+(S)*(x*x-rb*rb)/2
# S=wAw+kAk=d+2*wAw
assert s.factor(e_rel.subs(S,d+2*W)-(Ueff-Ueff.subs(x,rb))-W*(x*x-rb*rb))==0
print("homogeneous-background effective potential and relative-energy identity: verified")
Ebg=s.diff(rho,r,2)+s.diff(rho,r)/r-lam*rho*(rho*rho-rb*rb)
test=s.diff(r*s.diff(rho,r)*(rho-rb),r)-r*s.diff(rho,r)**2-r*lam*rho*(rho-rb)**2*(rho+rb)
assert s.simplify(s.expand(test-r*(rho-rb)*Ebg))==0
print("independent positive-amplitude multiplier identity: verified")
# Static surface balance from connection contraction on a cylindrical shell.
R,ps,jump=s.symbols("R ps jump", positive=True)
Gamma_rphiphi=-R
Sphiphi=ps/R**2
assert s.simplify(jump+Gamma_rphiphi*Sphiphi-(jump-ps/R))==0
print("cylindrical surface traction sign: verified")
print("No PDE solves, numerical integration, gravitational construction or certification.")
~~~

Observed output, with every assertion completed:

~~~text
SymPy 1.14.0 - exact localization identities
general cylindrical virial and off-shell force balance: verified
independent Cartesian transverse stress trace: verified
Derrick derivative: verified
homogeneous-background effective potential and relative-energy identity: verified
independent positive-amplitude multiplier identity: verified
cylindrical surface traction sign: verified
No PDE solves, numerical integration, gravitational construction or certification.
~~~

## 4. What these checks do not establish

No independently selected external reviewer has verified the package. Common-source transcription errors remain a possible limitation. The broad action leaves V and kinetic functions unspecified; this investigation does not prove nonexistence for every admissible choice, every field topology, every vacuum or nonstationary field. Pointwise positivity on a punctured/open domain is weaker than a uniform positive lower bound at its omitted endpoints.

The old exact matter solutions and Eq. (34) correction remain unchanged. No historical PASS, frozen Run 43 score, or guessed Gaussian-source data is used as localization proof. The new obstruction has its own assumptions and does not downgrade the established extended-matter existence verdict or upgrade the unestablished complete-spacetime verdict.
