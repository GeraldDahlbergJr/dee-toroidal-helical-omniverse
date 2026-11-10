# Stage F: non-inheriting stationary matter feasibility

2026-11-10 UTC. Exact parent `877cf5ae849279b9b2657a67122dc9959d34b233`.
One Codex agent; branch `research/stage-f-target-space-symmetry`.
The unchanged action, complete target classification and current conventions
are in `DEE_STAGE_F_TARGET_SPACE_SYMMETRIES.md`.

**Outcome B — precisely characterized residual class.** For distinct target
coupling eigenvalues, all nontrivial globally symmetry-generated stationary
motions fail the asymptotic fixed-vacuum necessary condition. For repeated
eigenvalues A=a I, real phase-plane rotations about the limiting vacuum evade
that obstruction. The metric enhancement at A=0 contributes no further
potential-preserving mixed motion. A surviving rotation is a candidate
representation, not an existing localized Einstein–DEE solution. The scalar
phase equation is massless at infinity and creates a further serious spectral
localization question. A general solution or a universal stationary no-go is
**NOT YET ESTABLISHED**.

Labels: **VERIFIED** means a proved statement under its stated hypotheses or
an identified exact local symbolic check, not independent peer review.
**CONDITIONAL** marks assumptions requiring separate validation.
**OBSTRUCTION IDENTIFIED** excludes a class only within those assumptions.
**NOT YET ESTABLISHED** never denotes evidence for existence.

## Stationary stress and equivariance

Let k be a spacetime Killing vector. With C_mu nu=G_IJ d_mu Phi^I d_nu Phi^J,

    T_mu nu=C_mu nu-g_mu nu[g^ab C_ab/2+V].

The condition L_k T=0 does not by itself assert L_k Phi=0. First restrict to
an equivariant map

    dPhi(k)=X(Phi),                                     (F8)

where X is a full-action target Killing vector on a target domain containing
the image. With local stationary time t and flow F_t of X, the complete local
representation is `Phi(t,x)=F_t(phi(x))`. The flow and map must stay in the
positive domain. The time derivative of a scalar gradient is
`L_k(d Phi^I)=partial_J X^I d Phi^J`. Therefore

    L_k C_mu nu=(L_X G)_IJ d_mu Phi^I d_nu Phi^J,
    L_k T_mu nu=P_mu nu-g_mu nu[g^ab P_ab/2+X(V)],
    P_mu nu=(L_X G)_IJ d_mu Phi^I d_nu Phi^J.             (F9)

Full-action symmetries make both terms zero. This proves stationary stress
kinematically without circularity, without Einstein equations, and without
assuming the scalar equations. It includes spatially varying rotational
profiles not representable as Stage E's Q=w t+m phi+q.

The scalar Euler derivative E_I transforms as a target covector under a
full symmetry, so E=0 at one time for the substituted profile is transported
to the other times. This covariance does not make E=0 automatic: the spatial
profile must solve all three reduced scalar equations, with the time derivatives
specified by X and `partial_t² Phi^I=partial_J X^I X^J` in adapted coordinates.
The Einstein equations additionally require the stationary metric to solve
`R_mu nu=8pi G_Newton(C_mu nu+V g_mu nu)` in every component. Bianchi/stress
conservation gives E_I d_nu Phi^I=0, which need not imply each E_I=0 when
the map has deficient rank. Neither (F8) nor (F9) licenses omitting equations.

For a metric Killing field broken by the potential, (F9) reduces to
`L_k T_mu nu=-g_mu nu X(V)` on such a map. Thus stationarity would require
X(V)=0 along its image, a weaker image condition than global preservation.
For this polynomial V a continuous trajectory staying on a fixed potential
level has constant rho: its permitted rho values at fixed V form a finite
set. This observation does not classify all possible image-restricted
motions of broken generators or all stationary stresses. The main theorem
below explicitly concerns potential-preserving target symmetries.

No proof is given that every map with L_k T=0 satisfies (F8) for a single
Killing field on the whole target, or even on an image neighborhood. Low-rank
maps, image-specific transformations and non-equivariant realizations cannot
be discarded using this analysis. This is an essential boundary on all
exclusion statements below.

## Observer energy and the vacuum fixed-point theorem

For any unit timelike observer u and orthonormal spatial frame e_a,

    epsilon_u=T(u,u)=|dPhi(u)|_G²/2
                         +sum_a |dPhi(e_a)|_G²/2+V.     (F10)

All terms are nonnegative in the authorized sector. A null Lorentz kinetic
contraction cannot cancel this energy (Stage E's null-gradient check remains
valid). In an asymptotic region where k is timelike, set
`F=sqrt(-g(k,k))`, u=k/F. Under (F8) the exact observer kinetic term is

    epsilon_u >= |X(Phi)|_G²/(2F²).                    (F11)

This requires no circularity and no assumption on spatial gradients. Timelike
k is needed only near infinity, not in the interior; ergoregions are not
excluded by (F11). Normalize the asymptotic inertial time so F->1 (a finite
positive constant would give the same obstruction with a rescaled bound).
A rotating time vector becoming spacelike or an unbounded lapse at infinity
would not meet these hypotheses.

For a spacelike slice, k=N n+beta with beta tangent, the exact alternative is

    Pi^I=dPhi^I(n)=(X^I-beta^i D_i Phi^I)/N,
    epsilon_n=(|Pi|_G²+h^ij G_IJ D_i Phi^I D_j Phi^J)/2+V,
    j_i=-G_IJ Pi^I D_i Phi^J.                          (F12)

A finite slice matter energy is E_m=integral epsilon_n dvol_h, not ADM mass.
To compare (F11) to this energy assume N->1, beta->0 in an asymptotically
orthonormal frame; then u and n have uniformly bounded relative boosts and
approach one another. Alternatively assume differentiated scalar falloff with
beta^i D_i Phi->0 and use (F12) directly. Spatial derivatives cannot rescue
nonzero X_infinity while retaining vanishing positive energy. For a bound
without derivative falloff, (F12) and Cauchy–Schwarz give

    |X|²=|N Pi+beta^i D_i Phi|²
       <=2 N² |Pi|²+2 |beta|_h² |DPhi|²,
    epsilon_n >= |X|²/[4 max(N²,|beta|_h²)]             (F13)

where the denominator is bounded above and nonzero on the normalized end.
This weaker bound suffices for energy divergence. The sharper limit from
(F11) gives the observer-frame coefficient. No cancellation of shifts or
spatial gradients has been presumed.

### Theorem F1: necessary fixed vacuum

**VERIFIED / OBSTRUCTION IDENTIFIED under explicit end hypotheses.** Suppose
an equivariant smooth scalar map into a positive target domain satisfies (F8)
for a smooth full-action generator X, on an asymptotically Euclidean end with
infinite proper volume. Assume V>=0, normalized stationary k timelike near
infinity, a regular asymptotic slice as above, and a limiting target value
`Phi->Phi_infinity=(v,Q_infinity)` uniformly on the end, with M(v)>0. Then:

1. Vanishing asymptotic orthonormal stress requires X(Phi_infinity)=0.
2. Finite E_m also requires X(Phi_infinity)=0, independently of pointwise
   stress decay, under these limiting-field and end hypotheses.
3. If X(Phi_infinity)!=0, epsilon_u has a positive lower bound outside a
   compact set and E_m diverges; differentiated asymptotically vacuum Einstein
   falloff is incompatible with that stress.

Proof: continuity makes |X(Phi)|_G² converge to a strictly positive number
if X_infinity!=0. Use (F11), F->1, and (F13) for the slice energy. The end has
infinite volume. Nonnegative other terms cannot cancel the bound. Conversely,
zero asymptotic stress makes the left side of (F11) tend to zero, and strict
G(v) positivity forces X_infinity=0. QED. These are necessary conditions only.
The observer identities are machine checked; the continuity, global end and
integral arguments are analytical.

If Q has no limiting vacuum value, the expression X(Phi_infinity) is not
available. Vanishing stress still requires |X(Phi)|->0 by (F11), but this is
not a theorem supplying a target limit. Translations remain obstructed
without a Q limit since their norm only depends on rho and M(v)>0. For a
nonzero plane rotation, its norm tends to zero only as Q approaches its
unique center, so in that case the condition itself forces the phase limit
in a real phase chart. Circle maps require convergence as maps in the target,
not as arbitrarily chosen discontinuous representatives. If vacuum limits
are angle-dependent, apply the pointwise condition at each limiting value;
finite energy additionally controls angular derivatives. No assertion that
finite energy alone enforces uniform pointwise decay is made.

Pointwise zero stress is weaker than finite energy. For the rotation residual
with frequency Omega!=0, (F11) near vacuum implies square integrability of
`Q-Q_c` on the end. A uniform power tail |Q-Q_c|=O(R^-p) gives integrable
temporal energy for p>3/2; a nonzero leading R^-p over a fixed solid-angle
belt fails it for p<=3/2. Spatial gradient energy alone only asks p>1/2 for
such differentiated tails. With rho-v=O(R^-s), V is quadratic near v, so
s>3/2 suffices for potential integrability and the differentiated amplitude
gradient integral. These sufficient power conditions are not solutions or
optimal universal necessary rates. Weighted/integrated decay, not just powers,
is the actual requirement. Proper and lapse-weighted matter energies are
comparable on this end; neither is silently identified with ADM mass.

### Theorem F2: classification of surviving generators

**VERIFIED / CONDITIONAL.** Under F1's hypotheses and a single full-action
one-parameter symmetry, the complete list is:

* If A has distinct eigenvalues, X=c^a d_Qa and X_infinity=0 forces c=0.
  Thus every such configuration inherits k; there is no nontrivial
  symmetry-generated non-inheriting motion in this class.
* If A=a I, including a=0, then X^Q=c+Omega J Q with
  J=[[0,-1],[1,0]], X^rho=0. If Omega=0, F1 again forces c=0.
  If Omega!=0, there is a unique fixed center `Q_c=J c/Omega` and

      X^Q=Omega J(Q-Q_c), Q_infinity=Q_c,
      |X|²=(1+a rho²)Omega² |Q-Q_c|²,
      Q(t,x)=Q_c+exp(Omega t J) q(x), rho=rho(x).       (F14)

  These are precisely the nonzero generators passing F1. They can vanish
  at vacuum while remaining nonzero in the interior. No mixed-rho generator
  is added when a=0 because V breaks that metric enhancement.
* If both phases are circle-valued and X is required to be a global symmetry
  of the entire torus target, only translations descend. F1 eliminates every
  nonzero such generator, for all A.

Proof: insert the complete full-action algebra from the symmetry report into
X_infinity=0. In two dimensions a nonzero skew generator is invertible; solve
for its center. Torus lattice descent removes it. QED. The conclusion for
circle phases is deliberately restricted to **global target symmetries**.
A small image confined to an embedded disk of the phase torus can use its
Euclidean chart and the locally defined rotational generator for A=a I.
For sufficiently small bounded q, (F14) stays in that disk for all t and
is a smooth torus map with stationary stress. This is a permissible kinematic
local-image realization, although rotation is not an action on arbitrary
torus fields. It is not excluded by the global-generator theorem and is not
an established solution. Local Killing classification and global descent
must not be conflated to claim a theorem for all torus stationary stresses.

F2 sets a generator to zero for the distinct case; it does not set all spatial
field gradients to zero, prove static geometry, or exclude every inheriting
noncircular stationary solution. Stage E's additional hypotheses are still
required for its Minkowski conclusion.

## Residual equations, regularity and Einstein constraints

Take the repeated case and define the notation
`Z=(Theta-Theta_c)+i(Psi-Psi_c)=exp(i Omega t) z(x)`, f=1+a rho².
This is notation for the unchanged two real fields, not a new complex-scalar
action or a polar redefinition of the physical phases. The exact equations are

    Box rho-lambda_V rho(rho²-v²)-a rho d_mu Q dot d^mu Q=0,
    nabla_mu(f d^mu Z)=0.                              (F15)

On any stationary metric in adapted coordinates, let H=sqrt(-det g). Then
the reduced complex equation, equivalent to both real phase equations, is

    (1/H) partial_i(H f g^ij partial_j z)
      +2 i Omega f g^ti partial_i z
      +(i Omega/H) partial_i(H f g^it) z
      -Omega² f g^tt z=0.                              (F16)

No circularity is assumed here; rho and f inherit k, and all metric coefficients
are stationary. The amplitude equation retains its time-dependent phase
kinetic contraction, which becomes time independent under (F14). The scalar
current integral used for static fields in Stage D now contains temporal
terms of opposite Lorentz sign and is no longer coercive.

For a stationary slicing with convention
`K_ij=(-partial_t h_ij+D_i beta_j+D_j beta_i)/(2N)`, all candidates must satisfy

    R[h]+K²-K_ij K^ij=16 pi G_Newton epsilon_n,
    D_j(K^ij-h^ij K)=8 pi G_Newton j^i,
    j_i=-G(Pi,D_i Phi),
    S_ij=G(D_i Phi,D_j Phi)+h_ij L_m,
    epsilon_n+S=2|Pi|²-2V.                             (F17)

The spatial Einstein evolution equations must hold as well. These constraints
and the full trace-reversed equations stated above are independent necessary
conditions, not a constructed constraint solution. Positive kinetic energy
requires f>0 on the entire profile, with a uniform positive lower bound on
compact smooth image ranges. If a<0 the profile cannot reach the finite-rho
boundary f=0; f(v)>0 alone says nothing about the interior. No sign of a is
otherwise needed for F1–F2.

Regularity is regularity of the original scalar maps in nonsingular spacetime
charts. There is no field-coordinate singularity at Q=Q_c: both Cartesian
phase fields are independent smooth scalars there. A rotation can therefore
have a fixed point although M(v)>0. This differs from an additive time drift,
whose vector never vanishes. Smooth bounded profiles with sufficiently fast
decay satisfy the kinematic positivity and finite-energy conditions, but
must still solve (F15)–(F17) and the remaining Einstein equations.

**Explicit boundary example (VERIFIED kinematics only):** on Minkowski space
choose rho=v, `q(x)=epsilon exp(-R²)(1,0)` and Omega!=0. This is smooth,
non-inheriting, has stationary stress and finite matter energy, fixes the
vacuum at infinity, and has no axis defect. For torus phases choose epsilon
small enough to stay in an embedded disk. It is not an on-shell counterexample
to any no-go theorem: `Delta exp(-R²)+Omega² exp(-R²)
=(4R²-6+Omega²)exp(-R²)` is nonzero, so its phase equations fail even at a=0;
its nonzero stress also fails the flat Einstein equation. This example shows
exactly why symmetry, localization and smoothness do not suffice for existence.
It prevents treating a purely kinematic obstruction as exhaustive.

### Filled axes and circular subsidiary calculation

No axis is needed in F1–F2. If axisymmetry is additionally imposed, a smooth
real phase doublet can have `Z=exp(i(Omega t+n phi)) F(s,z)` with integer n,
provided its Cartesian expression extends smoothly through the filled axis.
Here s is transverse proper distance. A sufficient local form is
`F=s^|n| f_0(s²,z)` with smooth f_0 (and smooth even rho); more general smooth
Cartesian modes obey corresponding Taylor restrictions. Then
`|partial_phi Q|²/W²=O(s^(2|n|-2))` and temporal energy is bounded; the local
volume s ds dphi dz gives finite axis energy for |n|>=1. At the axis Q=Q_c,
and dPhi(eta)=0. These are rotations of a Cartesian doublet about its center,
not nonzero degrees of the original independent circle phases. For n=0 the
usual smooth even scalar conditions apply. Stage E's Q=m phi singularity
and disk-degree obstruction remain intact; no zero of rho cures that ansatz.
Real Q imposes integer n by single-valuedness of this rotating profile, while
circle-valued phases in a small disk have no phase winding around their own
circle factors. Positive M(0)=I does not obstruct a smooth doublet zero.

For a **subsidiary circular** metric
`ds²=-N²dt²+gamma_AB dx^A dx^B+W²(dphi-omega dt)²`,
with real F=F(x^A), the preceding profile has

    Pi_Q=(Omega+omega n) J(Q-Q_c)/N,
    j_phi=-f (Omega+omega n)n |F|²/N,
    epsilon_temporal=f (Omega+omega n)² |F|²/(2N²),
    epsilon_angular=f n² |F|²/(2W²).

Because Jq dot q_A=0 for a real amplitude F, C_tA=C_phiA=0, so it meets
those algebraic circularity restrictions. For a general complex spatial z
with varying internal argument, these cross components need not vanish; the
circular Einstein equations then impose further restrictions. These are
kinematic compatibility examples, not solved geometries. Nonzero j_phi can
source rotation; Stage E's source-free rotation identity cannot simply set
omega=0. Even n=0 and vanishing axial momentum do not establish a general
noncircular staticity theorem. A nonzero Pi contributes positive active mass
in (F17), so Stage D's `epsilon+S=-2V` and its static mass proof cannot be
used until Pi and the relevant geometry have actually been removed.

## Massless asymptotic barrier and a limited additional theorem

With differentiated asymptotically inertial falloff, rho->v, f->f(v)>0,
(F16) has leading equation

    (Delta_flat+Omega²) z=0,                            (F18)

rather than a massive equation with exponentially decaying subthreshold
frequencies. The potential supplies radial mass
`m_rho²=V''(v)=2 lambda_V v²`, but zero phase masses. Small-amplitude massless
phase waves do not acquire a mass from their constant positive stiffness.
Ordinary nonzero radiative multipoles of (F18) have oscillatory 1/R tails,
whose temporal energy integral diverges. This is another **OBSTRUCTION
IDENTIFIED** for those tails, not a proof of their completeness on an arbitrary
coupled stationary end. A rapidly decaying profile cannot be justified merely
by replacing (F18) with a massive boson-star equation.

**Theorem F3 (VERIFIED, restricted reference problem).** On all of Minkowski
R^{1,3}, with A=0 and smooth time-harmonic phase doublet at real Omega!=0,
finite temporal energy implies z in L²(R³). The phase equations require
`(Delta+Omega²)z=0` distributionally. Fourier transformation gives
`(|p|²-Omega²) z_hat=0`. Thus the L² Fourier function is supported on the
sphere |p|=|Omega|, a measure-zero set, so z=0. At Omega=0 there is no
non-inheriting motion of this kind. This Fourier argument is analytical;
the verifier checks the local Helmholtz reduction, not functional analysis.
It does not prove nonexistence on a self-gravitating stationary spacetime,
for a nonconstant f, or even by itself in a curved exterior domain.

A theorem excluding square-integrable real-frequency modes of the *coupled*
operator (F16), with the actual interior geometry, coefficient decay and
possible ergoregions, would need its own proof and hypotheses. Leading
asymptotic balance, a 1/R generic tail, or a numerical sample is insufficient
to exclude embedded modes, a vanishing exterior with continuation issues,
or all allowed global interiors. No such spectral/unique-continuation theorem
is imported without proof. The residual sector is therefore **CONDITIONAL /
NOT YET ESTABLISHED** dynamically; it is not excluded by F1, not proved to
exist, and unresolved after this bounded stage.

## Scientific conclusion and limits

| Finding | Status |
|---|---|
| Complete local target metric and full-action symmetry classification | VERIFIED |
| No mixed-rho full-action generator, including A=0 | VERIFIED |
| Nonzero vacuum generator incompatible with normalized isolated end | VERIFIED / OBSTRUCTION IDENTIFIED under F1 |
| Distinct-eigenvalue equivariant non-inheriting motion | OBSTRUCTION IDENTIFIED under F1–F2 |
| Global torus symmetry-generated non-inheriting motion | OBSTRUCTION IDENTIFIED under F1–F2 |
| Repeated-eigenvalue rotations fixing a real vacuum or local image chart | CONDITIONAL residual representation |
| Smoothness, positive energy and integrated localization achievable kinematically | VERIFIED, not an on-shell solution |
| Flat free time-harmonic finite-energy phases | OBSTRUCTION IDENTIFIED under F3 |
| Coupled stationary residual solutions or their exclusion | NOT YET ESTABLISHED |
| Exhaustion of all stationary stress by target equivariance | NOT YET ESTABLISHED |
| Universal stationary no-go, stability and PRD readiness | NOT YET ESTABLISHED |

The strongest defensible result is F1–F2: a complete target-symmetry reduction
to a vacuum-centered rotational sector at A=a I, with an exact localization
obstruction for all generators lacking a vacuum fixed point. F3 explains why
the remaining motion is not automatically a viable boson-star mechanism.
Horizon boundaries, other ends and nonnormalized time generators lie outside
F1's setting; vacuum Kerr still warns against inferring global staticity
from zero matter flux (Stage E). No prior theorem is broadened silently.

The single most important unresolved question is whether the coupled stationary
Einstein–DEE equations admit or exclude a smooth finite-energy real-frequency
vacuum-centered rotational mode in the repeated-eigenvalue sector, under
precisely controlled asymptotic and global geometric hypotheses. An answer
would still not classify non-equivariant stationary stresses universally.

This work uses standard sigma-model geometry, Killing/Noether theory and
massless wave analysis. The model-specific classification and fixed-vacuum
obstruction are useful bounded results; novelty relative to the full literature
and independent peer validation are NOT YET ESTABLISHED. A conservative
Physical Review D assessment is **not ready**: no complete stationary
existence/nonexistence result for the residual sector, stability analysis,
independent mathematical review or exhaustive priority study is supplied.
No manuscript is drafted or submitted. Stop after Stage F.

## Reproduction and preservation record

Python 3, SymPy 1.14.0. The standalone verifier uses exact expressions, no
parameter sampling or imports of earlier verifiers. Named checks may contain
multiple tensor components; their count is not a count of independent global
theorems. Classification completeness uses the analytical curvature reduction,
not a black-box PDE solver. The verification is one investigator's independent
symbolic routes plus unchanged regressions, not an external reviewer team.

The following execution record is appended after actual runs. Logs live in
`/tmp`, outside the repository. Only the three authorized Stage F files are
added. Frozen-source supported mode is `--package-only`, as established in
Stages D/E; its historical full-tree addition guard is not a Stage F acceptance
criterion. The direct exact-parent tree/content comparison checks all previous
tracked files, including those outside the frozen package and all tracked
Run 43 artifacts. The unchanged provenance recovery gate remains INCOMPLETE
and physical tests remain NOT RUN; neither is represented as passed.

### Observed execution results, 2026-11-10 UTC

All eight commands below exited 0. Stage F reports **76 named exact local checks**, Stage E **61**, Stage D **27**, and Stage C **23**. The unchanged Stage A½ symbolic, A½ bounded and Stage B scripts each report their final PASS without a numbered check total; no invented count is assigned. A summary PASS line is not an additional check. Frozen provenance verifies **49 package hashes, 20 archive members and pinned source copies**. Earlier verifiers retain their original analytical scopes. The verifier was expanded from an earlier successful 68-check run to the final 76-check run; no FAIL occurred.

Actual complete stdout/stderr from the final runs follows (no output omitted):

Command: `python DEE_STAGE_F_VERIFIER.py`; exit code: **0**.

```text
SymPy 1.14.0
PASS: all target Christoffel components
PASS: direct Ricci scalar
PASS: curvature rational numerator
PASS: rank one curvature
PASS: opposite eigenvalues are not flat
PASS: repeated curvature
PASS: only flat parameter candidate from numerator coefficients
PASS: six unreduced Killing equations
PASS: translation x
PASS: translation y
PASS: rotation obstruction for distinct eigenvalues
PASS: repeated eigenvalue rotation
PASS: general affine phase generator obstruction
PASS: flat Killing T_rho
PASS: flat Killing T_x
PASS: flat Killing T_y
PASS: flat Killing L01
PASS: flat Killing L02
PASS: flat Killing L12
PASS: flat algebra closure T_rho,T_x
PASS: flat algebra closure T_rho,T_y
PASS: flat algebra closure T_rho,L01
PASS: flat algebra closure T_rho,L02
PASS: flat algebra closure T_rho,L12
PASS: flat algebra closure T_x,T_y
PASS: flat algebra closure T_x,L01
PASS: flat algebra closure T_x,L02
PASS: flat algebra closure T_x,L12
PASS: flat algebra closure T_y,L01
PASS: flat algebra closure T_y,L02
PASS: flat algebra closure T_y,L12
PASS: flat algebra closure L01,L02
PASS: flat algebra closure L01,L12
PASS: flat algebra closure L02,L12
PASS: phase translations commute
PASS: [R,T_x]=-T_y
PASS: [R,T_y]=T_x
PASS: potential action T_rho
PASS: potential action T_x
PASS: potential action T_y
PASS: potential action L01
PASS: potential action L02
PASS: potential action L12
PASS: phase affine potential preservation
PASS: mixed flat rotation potential breaking
PASS: off-shell arbitrary vector Noether identity
PASS: current contraction and sign T_x
PASS: current contraction and sign T_y
PASS: current contraction and sign R
PASS: current contraction and sign T_rho
PASS: current contraction and sign L01
PASS: rotation current components
PASS: stationary stress pullback kinetic identity
PASS: full stress Lie identity
PASS: stationary full action generator translations
PASS: stationary full action generator isotropic rotation
PASS: positive observer energy
PASS: ADM normal kinetic expansion
PASS: centered rotation norm
PASS: fixed vacuum centered rotation
PASS: translation vacuum norm
PASS: rotation center solves affine fixed point
PASS: Stage E additive normal velocity
PASS: Stage E M(0) stiffness
PASS: vacuum radial mass
PASS: vacuum phase masses
PASS: massless rotating asymptotic Helmholtz operator
PASS: Gaussian profile is not a free harmonic solution
PASS: noncircular stationary reduced phase equation
PASS: circular residual normal derivative
PASS: circular residual angular momentum
PASS: circular residual temporal energy
PASS: circular residual angular energy
PASS: circular residual mixed compatibility
PASS: finite Gaussian spatial L2 norm
PASS: rotation lattice discrepancy
PASS: Stage F 76 exact local checks; global proofs and limitations are in the reports.
```

Command: `python DEE_STAGE_A_HALF_VERIFIER.py`; exit code: **0**.

```text
SymPy 1.14.0
Einstein (0, 0) = 0
Einstein (0, 1) = 0
Einstein (0, 2) = 0
Einstein (0, 3) = 0
Einstein (1, 1) = 0
Einstein (1, 2) = 0
Einstein (1, 3) = 0
Einstein (2, 2) = 0
Einstein (2, 3) = 0
Einstein (3, 3) = 0
gauge 0 = 0
gauge 1 = 0
gauge 2 = 0
gauge 3 = 0
trace reversal involution = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
bar trace = -4*pi*G*(B - 4*epsilon)*(x**2 + y**2)
h = Matrix([[-2*pi*G*(B - 2*epsilon)*(x**2 + y**2), 0, 0, 0], [0, 2*pi*G*(B*x**2 + B*y**2 - 3*epsilon*x**2 - epsilon*y**2), -4*pi*G*epsilon*x*y, 0], [0, -4*pi*G*epsilon*x*y, 2*pi*G*(B*x**2 + B*y**2 - epsilon*x**2 - 3*epsilon*y**2), 0], [0, 0, 0, -2*pi*G*(B + 2*epsilon)*(x**2 + y**2)]])
direct linearized Einstein tensor = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
background amplitude = 0
background stress from action = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
both general phase equations from current density = [0, 0]
general amplitude equation from action = 0
all stress variations from action = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
conservation identity 0 = 0
conservation identity 1 = 0
conservation identity 2 = 0
conservation identity 3 = 0
frozen phase residuals = [0, 0]
frozen conservation residuals = [0, 0, 0, 0]
reported frozen amplitude residual = 0
frozen E_rho = -2*pi*G*rho_0*(B + 2*epsilon)*(x**2 + y**2)*(a*k_Theta**2 + b*k_Psi**2 + 2*c*k_Psi*k_Theta)
fixed-phase transverse amplitude reduction = 0
fixed-phase transverse phase equations = [0, 0]
fixed-phase transverse RHS = 2*pi*G*rho_0*(B + 2*epsilon)*(x**2 + y**2)*(a*k_Theta**2 + b*k_Psi**2 + 2*c*k_Psi*k_Theta)
PASS: all exact symbolic checks; no global existence or localization claim.
```

Command: `python DEE_STAGE_A_HALF_BOUNDED_VERIFIER.py`; exit code: **0**.

```text
SymPy 1.14.0
particular coefficients = [{alpha: -F/m**2, beta: -4*F/m**4}]
particular residual = [0]
I0 homogeneous residual = [0]
K0 homogeneous residual = [0]
I0 radial derivative = [0]
K0 radial derivative = [0]
general radial amplitude residual = [0]
I0 axis value = [0]
I0 axis derivative = [0]
K0 axis singularity = +oo
disk Dirichlet trace = [0]
disk Neumann derivative = [0]
disk Dirichlet radial residual = [0]
disk Neumann radial residual = [0]
disk zero-trace derivative obstruction expression = [0]
annular Dirichlet coefficient system = [0, 0]
annular inner Dirichlet trace = [0]
annular outer Dirichlet trace = [0]
annular Neumann coefficient system = [0, 0]
annular inner outward Neumann = [0]
annular outer outward Neumann = [0]
annular Dirichlet radial residual = [0]
annular Neumann radial residual = [0]
Green matching identity = [0]
positive-weight forcing integral primitive = [0]
integrated amplitude identity = [0]
disk Neumann mean formula = [0]
annular Neumann mean formula = [0]
inherited metric D^z on response = [0]
Theta general phase perturbation equation = [0]
Theta fixed-phase current residual = [0]
Psi general phase perturbation equation = [0]
Psi fixed-phase current residual = [0]
radial forcing agrees with inherited h_zz = [0]
full first-order fixed-phase radial amplitude residual = [0]
PASS: exact bounded response, boundary and both phase-current identities.
Analytical positivity/coercivity and matching proofs are in the companion documents.
```

Command: `python DEE_STAGE_B_VERIFIER.py`; exit code: **0**.

```text
SymPy 1.14.0
exterior Cauchy coefficient system = 0
exterior amplitude trace = 0
exterior normal derivative trace = 0
exterior original-action amplitude equation = 0
Abel Wronskian identity = 0
I0 small-radius normalization = 0
r K0 derivative normalization = 0
r I0 derivative K0 normalization = 0
Cauchy recovery of interior coefficients = 0
forcing sign from h_zz = 0
transverse cylindrical conversion = 0
direct extrinsic curvature coefficients = 0
Gauss normal Einstein constraint = 0
normal stress Einstein equation = 0
same-polynomial induced metric jump = 0
same-polynomial extrinsic curvature jump = 0
Theta current divergence = 0
Theta radial normal current = 0
Psi current divergence = 0
Psi radial normal current = 0
first-order scalar normal flux = 0
covariant conservation t = 0
covariant conservation r = 0
covariant conservation phi = 0
covariant conservation z = 0
background strong-energy combination = 0
PASS: Stage B exterior Cauchy matching, currents, junction geometry, constraint and conservation.
Finite collar only; analytical admissibility and obstruction scopes are in the companion files.
```

Command: `python DEE_STAGE_C_VERIFIER.py`; exit code: **0**.

```text
SymPy 1.14.0
PASS: rho variational equation vs covariant reduction
PASS: both phase variational equations vs current reduction
PASS: direct off-shell stress identity t
PASS: direct off-shell stress identity r
PASS: direct off-shell stress identity phi
PASS: direct off-shell stress identity z
PASS: normal traction formula
PASS: normal Einstein constraint: direct Ricci vs Gauss
PASS: normal Einstein constraint: explicit expression
PASS: phase taper mixed partials
PASS: invariant-gradient closure requires k prime zero
PASS: unit phase stiffness at zero amplitude
PASS: zero amplitude potential is positive
PASS: vacuum potential 1
PASS: vacuum amplitude equation 1
PASS: vacuum potential -1
PASS: vacuum amplitude equation -1
PASS: kinetic determinant polynomial
PASS: established axial amplitude branch
PASS: null phase contraction
PASS: null wave energy retains positive quadratic form
PASS: orthonormal energy sum of squares
PASS: radial vacuum fixed point
PASS: Stage C exact identities; analytical theorem scopes are in companion files.
```

Command: `python DEE_STAGE_D_VERIFIER.py`; exit code: **0**.

```text
SymPy 1.14.0
PASS: variational vs covariant rho
PASS: variational vs covariant Theta
PASS: variational vs covariant Psi
PASS: off-shell stress divergence t
PASS: off-shell stress divergence r
PASS: off-shell stress divergence phi
PASS: off-shell stress divergence z
PASS: Hamiltonian geometry G_nn=R3/2
PASS: lapse geometry R_tt=N D2N
PASS: momentum geometry 1
PASS: momentum geometry 2
PASS: momentum geometry 3
PASS: mixed Einstein geometry
PASS: static active mass eps+traceS=-2V
PASS: energy
PASS: phase integral identity
PASS: amplitude integration identity
PASS: positive-domain potential factor
PASS: mixed partial closure rho
PASS: mixed partial closure Theta
PASS: mixed partial closure Psi
PASS: axial gradient cannot be specified alone
PASS: vacuum V and Vprime 1
PASS: vacuum V and Vprime -1
PASS: nonvacuum zero amplitude
PASS: kinetic determinant
PASS: vacuum amplitude mass
PASS: Stage D 27 exact checks. Global hypotheses are not symbolic conclusions.
```

Command: `python DEE_STAGE_E_VERIFIER.py`; exit code: **0**.

```text
SymPy 1.14.0
PASS: inverse metric
PASS: determinant
PASS: kinetic contraction and all temporal-angular cross terms
PASS: metric variation T (0, 0)
PASS: metric variation T (0, 1)
PASS: metric variation T (0, 2)
PASS: metric variation T (0, 3)
PASS: metric variation T (1, 1)
PASS: metric variation T (1, 2)
PASS: metric variation T (1, 3)
PASS: metric variation T (2, 2)
PASS: metric variation T (2, 3)
PASS: metric variation T (3, 3)
PASS: action vs covariant equation rho
PASS: action vs covariant equation Theta
PASS: action vs covariant equation Psi
PASS: temporal current
PASS: angular current
PASS: meridional conserved current 0
PASS: meridional conserved current 1
PASS: off-shell stress divergence t
PASS: off-shell stress divergence r
PASS: off-shell stress divergence z
PASS: off-shell stress divergence phi
PASS: positive observer energy
PASS: matter angular momentum
PASS: stationary active mass
PASS: extrinsic curvature trace
PASS: extrinsic curvature norm
PASS: shift definition K (0, 0)
PASS: shift definition K (0, 1)
PASS: shift definition K (0, 2)
PASS: shift definition K (1, 1)
PASS: shift definition K (1, 2)
PASS: shift definition K (2, 2)
PASS: Hamiltonian Gauss identity
PASS: momentum Codazzi sign 0
PASS: momentum Codazzi sign 1
PASS: momentum Codazzi sign 2
PASS: frame dragging divergence
PASS: normal Ricci lapse identity
PASS: spatial Einstein geometry (0, 0)
PASS: spatial Einstein geometry (0, 1)
PASS: spatial Einstein geometry (0, 2)
PASS: spatial Einstein geometry (1, 1)
PASS: spatial Einstein geometry (1, 2)
PASS: spatial Einstein geometry (2, 2)
PASS: circular Ricci tA 1
PASS: circular Ricci phiA 1
PASS: circular Ricci tA 2
PASS: circular Ricci phiA 2
PASS: rotation Green identity
PASS: vacuum temporal energy coefficient
PASS: vacuum temporal stress isotropy
PASS: winding coefficient at zero amplitude
PASS: kinetic determinant
PASS: null-wave kinetic cancellation
PASS: null-wave energy retains positive quadratic
PASS: null-wave flux retains full stress
PASS: rotation source coercive when w=0
PASS: ADM angular momentum tail normalization
PASS: Stage E 61 exact local checks. Global proofs are in the reports.
```

Command: `python evidence/einstein_dee_sources/verify_evidence.py --package-only`; exit code: **0**.

```text
Git preservation NOT CHECKED in package-only mode
PASS: 49 package hashes, 20 archive members, pinned source copies
Recovery gate: INCOMPLETE. Physical tests: NOT RUN; approval required.
```

### Exact-parent preservation

A read-only Python comparison enumerated `git ls-tree -rz 877cf5ae849279b9b2657a67122dc9959d34b233`, read every baseline blob from the worktree (symlink targets as Git stores them), and recomputed SHA-1 of `blob <byte-length>\0<content>`. Actual output; exit **0**:

```json
{
  "baseline": "877cf5ae849279b9b2657a67122dc9959d34b233",
  "tracked_blobs_checked": 1117,
  "run43_named_blobs": 233,
  "mismatches": []
}
```

Thus **all 1,117 previous tracked blobs match byte for byte**, including **233 paths whose names contain run43 or run_43**, and all other Run 43 content regardless of naming because every baseline blob was checked. The original action, prior mathematical reports/verifiers and frozen evidence are included. This is independent of the provenance package-only mode.

`git diff --check` exited **0**, output empty. `git diff --exit-code 877cf5ae849279b9b2657a67122dc9959d34b233 --` before staging exited **0**, output empty (untracked additions are checked separately). `git status --short` exited **0** with exactly:

```text
?? DEE_STAGE_F_NONINHERITING_FEASIBILITY.md
?? DEE_STAGE_F_TARGET_SPACE_SYMMETRIES.md
?? DEE_STAGE_F_VERIFIER.py
```

The three authorized additions are staged and checked against the exact parent before commit; earlier-file comparison must stay empty. Remote commit and all three blob IDs are checked after push and reported in the delivery response. A commit cannot embed its own final hash. No Run 43 experiment or physical test is rerun. No Stage G is begun.
