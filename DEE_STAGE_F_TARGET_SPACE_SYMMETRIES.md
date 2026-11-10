# Stage F: target-space symmetries and currents

2026-11-10 UTC. Exact Stage E parent:
`877cf5ae849279b9b2657a67122dc9959d34b233`.
Branch: `research/stage-f-target-space-symmetry`. One Codex agent; bounded
analytical work. No numerical spacetime, theory change or Stage G.

**VERIFIED:** the complete local Killing algebra on every connected open
positive-kinetic target domain is classified below. The full-action algebra
is smaller in the flat case. **CONDITIONAL:** global symmetries depend on the
field domain and phase identifications. The companion feasibility report
proves a necessary localization condition and identifies an unresolved residual
class; it constructs no Einstein–matter solution.

## Preserved theory and scope

Authority: frozen `evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`,
Eqs. (3)–(9), (14)–(15). Read with both Stage D reports and its verifier,
both Stage E reports and its complete verifier, Stage A½ symbolic verification
and Stage C vacuum/no-go results as needed. Prior conclusions retain their
hypotheses. All older tracked files, original sources and Run 43 are unchanged.

Signature (-,+,+,+). The exact polynomial realization authorized here is

    S = integral sqrt(-g)[R/(16 pi G) - G_IJ d_mu Phi^I d^mu Phi^J/2 - V],
    Phi=(rho,Theta,Psi), G_target=diag(1,M(rho)),
    M=I+rho² A, A=A^T real, V=lambda_V(rho²-v²)²/4,
    G_Newton>0, lambda_V>0, v>0, M>0.

The symbol G_IJ is the target metric, not Newton's constant. No new field
identification or coupling is introduced. Target coordinate changes used in
proofs do not replace the physical fields. Work on an open connected component
of `D={rho: I+rho² A>0}` times a phase chart. A nonnegative amplitude may have
a boundary at rho=0; Killing equations are first solved in its interior.
Degenerate endpoints where M loses positivity are excluded. Potential vacuum
rho=v must belong to D. For a signed scalar -v is another vacuum; the same
classification holds there. Discrete symmetries are not continuous generators.

## Killing equations and curvature reduction

An orthogonal change of phase coordinates diagonalizes A, with real eigenvalues
alpha,beta. Locally put Q=(x,y) in that basis,

    f=1+alpha rho²>0, h=1+beta rho²>0,
    ds_T²=d rho²+f dx²+h dy², X=U d_rho+P d_x+H d_y.

Diagonalization on a torus transforms its lattice too; it need not be a torus
automorphism in the original periods. This local calculation does not assume
that it is one. Nonzero Christoffels are

    Gamma^rho_xx=-alpha rho, Gamma^rho_yy=-beta rho,
    Gamma^x_rho x=Gamma^x_x rho=alpha rho/f,
    Gamma^y_rho y=Gamma^y_y rho=beta rho/h.

The six equations L_X G=0 are exactly

    U_rho=0,
    U_x+f P_rho=0, U_y+h H_rho=0,
    alpha rho U+f P_x=0, beta rho U+h H_y=0,
    f P_y+h H_x=0.                                      (F1)

No assumed polynomial form for X is used. To eliminate possible mixed fields,
use the sectional curvatures in the orthonormal target frame:

    K_rho x=-alpha/f², K_rho y=-beta/h²,
    K_xy=-alpha beta rho²/(fh),
    R_T=-2[alpha/f²+beta/h²+alpha beta rho²/(fh)].        (F2)

These follow from the connection; the verifier independently computes Ricci
and the scalar. A Killing field preserves R_T, hence U R_T'=0.

Here is a completeness argument valid also for negative or zero eigenvalues.
The numerator of `-R_T f²h²/2` is

    n(rho)=alpha+beta+5 alpha beta rho²
                +2 alpha beta(alpha+beta)rho⁴+alpha² beta² rho⁶.

If both eigenvalues are nonzero, the rational function R_T tends to zero as
formal rho tends to infinity. Were it constant on any open positive domain,
a polynomial identity would make it constant as a rational function, so that
constant would be zero. But n identically zero would require alpha+beta=0
and alpha beta=0, contradicting both being nonzero. This use of rational
identities does not extend the physical domain through a pole. With just one
nonzero eigenvalue, R_T=-2 alpha/(1+alpha rho²)² (or interchange), which is
nonconstant. Therefore for every A!=0, R_T' is nonzero on a dense open set
of each allowed interval. Smoothness forces U=0 everywhere, including rho=0
and isolated critical-curvature slices. Restricting fields to one rho slice
would not classify isometries of the three-dimensional target.

Now (F1) gives P_rho=H_rho=P_x=H_y=0. Locally P=P(y), H=H(x), and
`f P'(y)+h H'(x)=0`. This identity is polynomial in rho²; its constant
coefficient gives P'+H'=0. Independent x,y then make each derivative constant.
Thus

    P=c-Omega y, H=d+Omega x,
    (alpha-beta)Omega=0.                               (F3)

This proves completeness without guessing that translations and rotations
exhaust the algebra. In particular rank-one A, opposite-sign eigenvalues,
alpha+beta=0, and any negative repeated eigenvalue inside its allowed interval
produce no hidden mixed generators. For a real symmetric 2x2 matrix, repeated
eigenvalues imply A=a I in the original phase basis.

## Complete classification, including flat enhancement

| Coupling | Complete local metric Killing algebra | Full-action algebra |
|---|---|---|
| alpha!=beta (necessarily A!=0) | R²: T_x=d_x, T_y=d_y | same |
| A=a I, a!=0 | e(2): T_x,T_y,R=x d_y-y d_x | same |
| A=0 | e(3): three translations and three rotations in (rho,x,y) | e(2): T_x,T_y,R |

All statements in this table are **VERIFIED** local statements on an open
positive domain. The distinct case includes one zero eigenvalue. The repeated
case includes positive and negative a, with `1+a rho²>0`; asymptotic vacuum
additionally requires `1+a v²>0`. No nonzero A is locally flat or of constant
sectional curvature: nonconstant scalar curvature already rules out both.
There is no other exceptional enhancement in the stated polynomial family.
At A=0 the connection vanishes and (F1) becomes the Euclidean Killing system.
Differentiating it and permuting indices yields `partial_i partial_j X_k=0`:
add the differentiated (i,k) and (j,k) equations and subtract the differentiated
(i,j) equation. Hence

    X=c^i d_i+B^i_j Phi^j d_i, B^T=-B,
    Phi=(rho,x,y),                                      (F4)

all six generators, including T_rho and rotations mixing rho with phases.
They are genuine local metric isometries, missed by a phase-plane-only guess.

With `[X,Y]^i=X^j partial_j Y^i-Y^j partial_j X^i`,

    [T_x,T_y]=0, [R,T_x]=-T_y, [R,T_y]=T_x.

For flat e(3), `T_i=d_i`, `L_ij=Phi^i d_j-Phi^j d_i`:

    [T_k,L_ij]=delta_ki T_j-delta_kj T_i,
    [L_ij,L_kl]=delta_jk L_il-delta_ik L_jl
                       -delta_jl L_ik+delta_il L_jk.

Use L_ji=-L_ij. All 15 pairs of the six flat generators and the phase brackets
are checked for closure in the verifier.

## Potential preservation and global symmetries

A symmetry of the matter action must additionally obey

    X(V)=lambda_V rho(rho²-v²) U=0                       (F5)

identically on an open target domain, not just at a vacuum or critical point.
Since lambda_V,v>0, (F5) implies U=0 by continuity. For A!=0 this already
holds; for A=0 it removes T_rho and both mixed rotations from (F4). Cancellations
between a radial translation and mixed rotations cannot make
`c^rho+B^rho_x x+B^rho_y y` vanish for all x,y unless all three coefficients
vanish. Thus the full-action algebra is precisely translations plus those
phase rotations commuting with A: `X^rho=0, X^Q=c+BQ`,
`B^T=-B, [A,B]=0`. A metric Killing vector broken by V is not a conserved
Noether symmetry even though X(V) happens to vanish at rho=v.

For real Q in R² and the full allowed rho component, phase translations and
these rotations are globally defined complete flows. In the flat signed
rho target R³ the metric e(3) flows are global, but only e(2) preserves V.
If rho>=0, radial translations and mixed rotations generally do not preserve
the half-space, so they are only local metric fields; the full-action e(2)
still acts globally. Arbitrary additional target boundaries or restrictions
can reduce the global algebra; none is assumed silently.

For circle-valued phases with a rank-two lattice Lambda, a vector field on
the cover must satisfy `X(rho,Q+ell)=X(rho,Q)` for every ell in Lambda.
The affine generator obeys this iff B ell=0. A full-rank lattice forces B=0.
Hence the identity component of the global full-action group is just the
phase torus of translations, for every A. A finite lattice-compatible rotation
may survive discretely; it supplies no continuous non-inheriting generator.
At A=0 mixed Euclidean rotations similarly fail periodicity; for a signed
unrestricted rho the global metric identity component is radial translation
times phase torus, with radial translation again broken by V. These statements
also distinguish a flat torus's translations from the rotations of its cover.
If only one phase is periodic, a nonzero two-dimensional skew B still cannot
annihilate its lattice vector; the full-action rotation is lost. A map whose
image lies in a small phase chart may nevertheless have a local rotational
realization: this is discussed separately in the feasibility report. It is
not a global action symmetry on the phase torus.

## Noether identities with explicit signs

Let the covariant Euler derivative be

    E_I = G_IJ[Box Phi^J+Gamma^J_KL d_mu Phi^K d^mu Phi^L]
                 -partial_I V,
    delta S_m = integral sqrt(-g) E_I delta Phi^I + boundary.

These are lower target indices. In particular

    E_rho=Box rho-V'-d_mu Q^T M' d^mu Q/2,
    E_Q=nabla_mu(M d^mu Q).

Choose the current sign to agree with the original phase currents:

    J_X^mu=G_IJ X^I d^mu Phi^J.                         (F6)

The canonical current for the convention delta Phi=+epsilon X is -J_X;
this overall sign changes no conservation law. A local parameter gives
`delta L_m=-J_X^mu partial_mu epsilon
 -epsilon[(L_X G)_IJ d_mu Phi^I d^mu Phi^J/2+X(V)]`.
Combining with the Euler variation, or directly differentiating (F6), gives

    nabla_mu J_X^mu = X^I E_I + X(V)
                      +(L_X G)_IJ d_mu Phi^I d^mu Phi^J/2.  (F7)

Thus on shell a metric Killing field has divergence **+X(V)** in this current
convention; only a potential-preserving Killing field is conserved. This sign
also follows from flat X=T_rho: Box rho=V' gives div(d rho)=V'. The verifier
checks (F7) for arbitrary vector fields using unrestricted first and second
jets in a spacetime normal frame. Because this is tensorial, those local
identities give the covariant identity, without assuming Minkowski spacetime.

For every full-action generator the explicit currents are

    J_c^mu=c^T M d^mu Q,
    J_R^mu=(BQ)^T M d^mu Q   (or B(Q-Q_c) about a chosen center),
    J_{c+BQ}^mu=(c+BQ)^T M d^mu Q.

The first line contains the two original directed phase currents (7) of the
source. For A=a I and unit R, the additional current is
`(1+a rho²)(x d^mu y-y d^mu x)`. It is an internal rotation current, with
the field values as weights, not an angular momentum associated to spacetime
rotation and not a third independent translational transport channel.
It is globally conserved only when the target symmetry and map are global
(or consistently on a suitable local image domain). Charges require further
integrability and boundary-flux hypotheses; local divergence alone does not
guarantee a finite conserved integrated charge.

For completeness the flat, potential-broken metric generators have

    J_T_rho^mu=d^mu rho,             div J=V',
    J_L_rho x^mu=rho d^mu x-x d^mu rho, div J=-x V',
    J_L_rho y^mu=rho d^mu y-y d^mu rho, div J=-y V'

on shell. Phase generators remain conserved. This lists every generator of
the exceptional metric algebra and its potential status; linearity covers
all combinations.

## Verification boundaries

**VERIFIED by exact symbolic local identities:** connection, Ricci scalar,
Killing equations and listed generators, algebra closure, potential variation,
Noether identity and contractions. **VERIFIED by analytical proof:** curvature
reduction and completeness, and lattice descent under the stated target model.
No symbolic run certifies global topology or existence. **CONDITIONAL:** applying
any global symmetry to a selected restricted field domain. **NOT YET ESTABLISHED:**
exhaustion of stationary-stress scalar maps by this symmetry-generated subclass.
The companion report records the actual test outputs and preservation checks.
