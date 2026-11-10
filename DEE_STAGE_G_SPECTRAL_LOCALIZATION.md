# Stage G: spectral localization and coupled limitations

2026-11-10 UTC. Parent `767783984a8502d892cf52048cd43eadde79e552`,
tree `0da0d813d76108bc880d2e34c47f82b261026cd1`. One agent, three new files.
Notation and complete equations (G1)–(G10) are in
`DEE_STAGE_G_COUPLED_STATIONARY_EQUATIONS.md`.

**Outcome B — CONDITIONAL obstruction**, with **Outcome C** for the general
stationary sector. A real-frequency finite-energy rotation is excluded in the
explicit spherical-exterior and globally elliptic class of Theorem G2 below.
Gravitational localization in an arbitrary stationary Einstein–DEE geometry
is **NOT YET ESTABLISHED**, and neither is its general exclusion. The attractive
Coulomb term does not create a mass gap. No localized physical solution is
constructed. The following are analytical results, not independent peer review.

## Energy, falloff and regularity: necessary conditions

Assume a smooth connected complete slice, one normalized asymptotically
Euclidean end, no physical inner boundaries, singular shells or horizons,
N>0, smooth positive h, f>0 on the entire map, rho->v, z->0. Interior
compact sets then have positive lower bounds on N and f; on the end f->f_v
=1+a v²>0. Smoothness, positive energy, finite ADM mass, vanishing asymptotic
stress and finite proper matter energy are separate requirements.

**Theorem G1 (VERIFIED necessary integrability).** If E_m in (G6) is finite,
Omega!=0, N and |beta|_h are bounded on the end, and h is asymptotically
Euclidean, then z and Dz belong to L² on the end (hence globally if the
interior is compact and smooth). Indeed

    |Omega z|²=|N p+beta.Dz|²
        <=(N²+|beta|_h²)(|p|²+|Dz|_h²),
    epsilon_Z >= f Omega² |z|²/[2(N²+|beta|_h²)].          (G11)

The denominator is bounded above and f bounded below there; also
`epsilon_Z>=f|Dz|²/2`. This proof retains shift cancellation exactly; no
separate assumed pointwise derivative falloff is needed. Finite ADM mass by
itself is not a substitute for finite E_m. Conversely z in H¹ is sufficient
for the mode energy on a globally regular compact interior plus this bounded
end, but does not solve the metric or amplitude equations. The amplitude
potential near v is `lambda_V v² (rho-v)²+O((rho-v)³)`.

At an interior point, rho,z must be smooth in regular spatial charts, with
f never reaching zero (especially important when a<0). At an ergosurface
N>0 and h>0 can stay smooth; the original wave equation and (G6) remain
regular, while elliptic continuation of the reduced equation loses its
usual justification. At infinity, smoothness and integrability alone do
not supply uniform differentiated coefficient expansions.

If axisymmetry is additionally imposed and z=F(s,y) exp(i n phi) in a
regular filled-axis chart, n is an integer and smoothness requires the full
Cartesian complex field to extend smoothly. For a pure azimuthal Fourier
mode a sufficient local form is F=s^|n| times a smooth function of s²,y.
The normal energy includes f n²|F|²/(2s²); unlike Stage E's additive phase
winding, this Cartesian doublet can vanish at the axis. This supplies no
existence theorem. rho and invariant metric functions must have regular
even expansions and the metric must satisfy elementary flatness. No axis
is treated as a physical boundary. All fields in the primary problem are
smooth and there are no reflecting surfaces.

## Controlled end expansions and what is actually assumed

Write O_k(R^-p) for bounds on spatial derivatives through order k, with
each derivative gaining one inverse power, in specified Cartesian end
coordinates. Mere limits N->1,h->delta,beta->0 and finite ADM mass do not
imply a monopole expansion with an O(R^-2) differentiated remainder.
Stationary Einstein equations, asymptotic gauge, matter decay and absence
of an asymptotic boost must justify such expansions in any proposed solution.
No deduction of all those facts from finite E_m is claimed here.

For a useful controlled subclass assume a stationary rest end with

    N=1-m/R+O_2(R^-2),
    h_ij=(1+2m/R)delta_ij+O_2(R^-2),
    beta^i=O_2(R^-2), rho-v=O_2(R^-s), s>3/2,
    m=G M_ADM.                                         (G12)

The h expansion explicitly fixes the isotropic monopole gauge. The ADM
surface formula gives M_ADM=m/G for this h expansion. The identical lapse
coefficient is an additional stationary mass-equality/end hypothesis, as is
appropriate to a controlled rest-frame Einstein end; it is not inferred
from the h expansion alone. For example a Schwarzschild vacuum end in
isotropic coordinates has N=(1-m/(2R))/(1+m/(2R)),
h=(1+m/(2R))^4 delta and satisfies (G12). This example checks coefficients;
it is not an allowed horizonless interior or a replacement DEE solution.
In an areal-radius gauge h has anisotropic 1/R coefficients. Arbitrary
stationary coordinates cannot inherit (G12)'s coefficient-by-coefficient
formula unchanged. Finite ADM mass is not by itself sufficient for (G12).

Under (G12), with eta=rho-v,

    f=f_v+2a v eta+a eta²,
    w=N sqrt(h)=1+2m/R+O_2(R^-2),
    P^ij=(1-2m/R)delta^ij+O_2(R^-2),
    w f P^ij=f_v delta^ij+2a v eta delta^ij
                         +O_2(R^-2)+O_2(R^-s-1),
    w f/N²=f_v(1+4m/R)+2a v eta
                         +O_2(R^-2)+O_2(R^-s-1),
    w f beta^i/N²= f_v beta^i+O_1(R^-3)+O_1(R^-s-2).

Errors with f_v dimensional factors are understood coefficientwise; constants
are fixed. In (G1) multiplied by w/f_v the order 1/R correction to the
second-derivative density cancels. Thus the coordinate-normalized equation is

    Delta z+Omega²(1+4m/R)z
       + partial_i(B^ij partial_j z)+Omega² C z
       + (2a v/f_v)[partial_i(eta partial^i z)+Omega² eta z]
       + 2i Omega (beta^i+O_1(R^-3)+O_1(R^-s-2))partial_i z
       + i Omega partial_i(beta^i+O_1(R^-3)+O_1(R^-s-2))z=0, (G13)

where B^ij,C=O_1(R^-2)+O_1(R^-s-1)+O_1(R^-2s), with higher regularity
inherited when available. This is a divergence-form perturbation, not just
a multiplicative potential. Multiplying by a positive function has changed
the natural formal-symmetry measure, so it does not assert symmetry in flat
L². The leading shift derivative coefficient is O(R^-2), its divergence
O(R^-3). A rest-frame spin tail, when separately justified, is

    beta=-2G (J cross x)/R³+O_1(R^-3),
    2i Omega beta.D=-4i G Omega J.L/R³,
    L=(J cross x).D/J for J!=0.

For J along z and an azimuthal exp(i n phi), the displayed term becomes
`+4G Omega J n/R³`; it changes sign with Omega J n. This spin normalization
is a conditional stationary spin gauge, not a consequence of finite mass.
No omitted shift term is declared zero in the general problem.

A particularly transparent exact specialization has beta=0,
h=H(R)² delta, N=N(R), f=f(R). Put C0=NHf and B=log C0. Then

    Delta z+B' partial_R z+Omega² H²/N² z=0,
    u=sqrt(C0) z,
    Delta u+[Omega² H²/N²-(Delta B)/2-|grad B|²/4]u=0.    (G14)

Here (G12) means H=1+m/R+O_2(R^-2), so
`H²/N²=1+4m/R+O_2(R^-2)`. If f-f_v=O_2(R^-s),
B=log f_v+O_2(R^-min(2,s)). Derivative potential terms are
O(R^-min(4,s+2)). The transformed Schrödinger notation at a fixed frequency is

    [-Delta+V_eff(Omega)]u=Omega² u,
    V_eff=Omega²(1-H²/N²)+(Delta B)/2+|grad B|²/4
          =-4m Omega²/R+O(R^-2).                        (G15)

V_eff is NOT the original scalar V(rho). For positive ADM mass the leading
term is attractive, for negative mass repulsive, for zero mass absent.
The subleading curvature and amplitude terms have no universal sign; shift
terms are differential and sign-indefinite. Positive-kinetic matter with
V>=0 satisfies the dominant energy condition, so appropriate complete
asymptotically flat initial-data positive-mass hypotheses imply M_ADM>=0;
this is a global input, not a pointwise consequence of N alone. No bound
state follows from attraction: Omega²>0 lies in the massless continuum and
the potential itself depends on Omega. An attractive Coulomb Hamiltonian's
negative-energy levels do not establish these positive-frequency modes.

The exact amplitude perturbation equation, with eta=rho-v and P as before, is

    [w^-1 partial_i(w P^ij partial_j)-mu²]eta
       =lambda_V(3v eta²+eta³)+a(v+eta)(Xz-|p|²),
    mu²=2lambda_V v².                                  (G16)

Consequently there is a massive homogeneous amplitude scale, but a forced
phase bilinear source. An exponential Yukawa tail cannot be asserted when
this source has only algebraic decay. The coupling changes f and derivative
coefficients in (G13)–(G15); it does not generate a nonzero asymptotic phase
mass while rho->v and f_v>0. Bounds on eta derivatives in (G12) are assumptions,
not consequences of L²(z) alone. Oscillatory or poorly controlled derivatives
of coefficients are a genuine gap in an automatic spectral argument.

## A controlled exclusion proved without a short-range Coulomb assumption

**Theorem G2 (CONDITIONAL / OBSTRUCTION IDENTIFIED).** Let the full unchanged
Einstein–DEE equations hold in the smooth horizonless setting of G1 with
Omega real and nonzero. In addition assume:

1. Outside a compact set there is an exactly spherical coefficient end:
   beta=0, h=H(R)² delta, N=N(R), rho=rho(R), with H,N positive,
   N=1-m/R+O_2(R^-2), H=1+m/R+O_2(R^-2),
   f=f_v+O_2(R^-s), s>3/2. No spherical symmetry of z is assumed.
2. On the entire connected slice |beta|_h<N, with smooth stationary
   coefficients and f>0. Ellipticity is locally uniform on compact sets.
3. E_m<infinity and no physical boundary; z is smooth globally.

Then z=0. These extra geometric hypotheses have NOT been derived for all
stationary Einstein–DEE configurations. In particular exact spherical
coefficients on an exterior are restrictive in a rotating problem.

Proof: by G1, z in H¹ intersect L². Expand on standard spherical harmonics
on the exterior: z=sum z_lm(R)Y_lm. Smooth radial coefficients allow each
projection to satisfy its scalar ODE. With B=log(NHf), set
`y_lm=R sqrt(NHf) z_lm`. Direct differentiation of (G14) gives

    y''+[Omega² H²/N²-l(l+1)/R²
                    -B''/2-B'/R-(B')²/4]y=0,
       = y''+[k²+c/R+r_l(R)]y=0,
    k=|Omega|>0, c=4m k², r_l=O(R^-2) in L¹(dR).         (G17)

Also L² of z makes each y_lm in L²(dR), since NHf tends to a positive
constant and h volume is comparable to Euclidean volume. This step does
not exchange a large-l limit with an asymptotic expansion: each fixed l
is treated separately.

Here is an explicit long-range ODE argument. Let
`theta=kR+(c/(2k))log R`, w_±=exp(±i theta). Substitution into
`d²/dR²+k²+c/R` gives residuals
`[-c²/(4k² R²) ∓ i c/(2k R²)] w_±`, which are L¹ errors.
The 2x2 fundamental comparison matrix with columns (w_±,w_±') has bounded
matrix and inverse for large R; its determinant is -2i theta', bounded
away from zero. Writing (y,y')=W a yields a'=E(R)a with E in L¹, including
r_l. The integral Gronwall bound gives a finite limiting vector a_infinity;
the fundamental matrix for a'=E a has a finite invertible limit (its inverse
also solves an integrable matrix equation). Thus a_infinity=0 iff y=0.
For a nonzero solution,

    y=A exp(i theta)+B0 exp(-i theta)+o(1),
    |A|²+|B0|²>0.

The integral of the cross oscillation exp(2i theta) is o(T) on [R0,T]
(in fact bounded up to a finite correction, by integration by parts using
theta'->k and theta''=O(R^-2)). The o(1) remainder has o(T) integrated
contribution. Therefore
`T^-1 integral_R0^T |y|² -> |A|²+|B0|²>0`, contradicting L².
Every harmonic vanishes, so z=0 on the exterior.

To continue inward, (G1) has a smooth real symmetric positive definite
principal matrix fP, with locally bounded complex first/zero-order terms.
The standard local elliptic unique-continuation estimate applies: in a
compact coordinate ball the real and imaginary parts obey a common real
uniformly elliptic scalar principal operator with residual bounded by
C(|grad u|+|u|) for the two-component vector u=(Re z,Im z). The usual
Carleman inequality can be summed over these two components and absorbs
these lower-order couplings. Smooth coefficients exceed the Lipschitz
principal/bounded lower-order regularity requirement. Weak continuation
from a nonempty open zero set follows (the stronger vanishing-to-infinite-
order result is unnecessary). Connectedness propagates the zero set through
overlapping balls. This is the standard elliptic unique-continuation input,
not a new theorem proved by SymPy; see N. Aronszajn, J. Math. Pures Appl.
36 (1957), 235–249, unique continuation for second-order elliptic inequalities.
All required local ellipticity/regularity and diagonal-system conditions have
been checked above. Hence z=0 globally. QED.

The theorem proves an obstruction on any background with these properties,
so also on coupled solutions; it does not need to solve the interior Einstein
equations first. Mass m of either sign is allowed in the ODE proof. If an
exterior-vanishing solution admits a separate continuation theorem through
an interior ergoregion, that theorem could replace hypothesis 2, but none
is established here. A surface R=R0 in this proof is an artificial matching
surface, not an imposed reflecting inner boundary. No physical boundary is
introduced. The case Omega=0 is not covered by the positive-frequency ODE
argument and is not non-inheriting rotation.

## General spectral question: exact limits of the argument

For a general stationary end the long-range 1/R frequency term is not L¹
after radial reduction and cannot be put in a short-range Rellich argument.
G17 handles it by the explicit logarithmic Coulomb phase. Standard flat
Rellich lemmas concern exact exterior Helmholtz equations with their decay
or radiation conditions; their hypotheses fail for (G13) as written.
No general absence-of-positive-eigenvalues theorem is invoked here without
coefficient, derivative, domain, tensor-metric and long-range assumptions.

The spherical separation proof does not handle arbitrary anisotropic
second-order O(R^-2) metric perturbations and first-order shift couplings.
Such controlled classes may admit stronger long-range scattering theorems,
but establishing the needed operator domain and precise hypotheses is a
remaining task, not a counterexample. Without coefficient derivative control,
leading Helmholtz asymptotics alone do not rule out embedded eigenfunctions:
oscillatory long-range potentials in other spectral problems illustrate why
one must prove, rather than guess, absence of embedded spectrum. No such
potential is claimed to be realizable by the coupled Einstein–DEE system.

For the general operator pencil, the asymptotic continuum threshold remains
Omega²=0. Finite E_m forces positive-frequency L² profiles; neither the
Coulomb correction nor an attractive compact interior potential alone proves
or excludes an embedded mode. Interior trapping changes scattering and
resonances, but is not itself evidence of a real L² eigenmode. In the globally
elliptic controlled class, exterior exclusion plus unique continuation makes
any compact interior irrelevant to this conclusion. With an ergoregion, the
stationary equation can be of mixed type; this continuation step cannot be
silently imported. Frame dragging contributes real-frequency differential
terms, may destroy global ellipticity, and in a controlled rest end has
shorter range than the mass monopole. Its smallness at infinity alone does
not settle the interior or a general anisotropic spectral theorem.

## Coupled identities and boundary terms

All following identities require the stated integrals to exist. Write
F^i=f(P^ij D_j z+i Omega beta^i z/N²). From (G3), multiplication by z*
and taking the real part gives exactly

    integral N f (Xz-|p|²)dvol_h
       =lim integral_boundary N Re(z*F^i)s_i dA.         (G18)

The imaginary shift boundary piece here is zero, but interior cross terms
are retained inside |p|². For H¹ modes on a bounded normalized end, use
cutoffs chi_R with |Dchi_R|<=C/R: the boundary replacement term is bounded
by C/R times tail L²(z)L²(Dz), plus a bounded shift tail, and tends to zero.
This establishes the zero-flux identity by cutoff, rather than assuming
a pointwise spherical flux limit for an arbitrary H¹ function. It says
`integral Nf Xz=integral Nf |p|²`, not that either side vanishes.

Multiplying (G2) by eta=rho-v gives

    integral N[P^ij rho_i rho_j+eta V'
                   +a rho eta(Xz-|p|²)]dvol_h
      =lim integral_boundary N eta P^ij rho_j s_i dA.   (G19)

Under (G12) this boundary is O(R^(1-2s))->0; regular interior charts have
no physical boundary. Even for rho>=0, where eta V'>=0, the phase source
has no fixed pointwise sign; when a rho eta/f varies, (G18) does not cancel
it. With an ergoregion even the first gradient term is indefinite.

Integration of (G10), using D_i beta^i=NK, yields

    lim integral_boundary D_s N dA
      =integral N[K_ij K^ij-K²+8 pi G(U-V)]dvol_h
                         +lim integral_boundary K beta.s dA. (G20)

The momentum constraint contracted with beta gives

    integral N(K_ij K^ij-K²)dvol_h
       =lim integral_boundary beta^i(K^j_i-delta^j_i K)s_j dA
                              -8 pi G integral beta^i j_i dvol_h. (G21)

For beta=O_1(R^-2), K=O(R^-3), these shift surface terms vanish as
O(R^-3). Smooth filled-axis excision terms vanish in regular Cartesian
geometry. No unmentioned inner boundaries are omitted. If lapse flux
equals 4 pi G M_ADM (the stationary end mass equality), then

    M_ADM=integral [N(epsilon+S)-2 beta.j]dvol_h
          =2 Omega Q_rot-2 integral N V dvol_h,
    Q_rot=integral f Im(z* p)dvol_h.                     (G22)

The second equality uses Np+beta.Dz=i Omega z: amplitude q contributions
cancel between NU and beta.j. Q_rot is the chosen positive normal internal
rotation charge; its sign convention is explicit. It is not spacetime
angular momentum. Under G1 and bounded coefficients it is finite. The
local current is f Im(Z* d^mu Z); charge conservation additionally needs
vanishing flux. No definite sign of Omega Q_rot is assumed for arbitrary
shift. Even beta=0 gives the positive contribution
`2 Omega² integral f |z|²/N`; Stage D's nonpositive static mass identity
therefore fails at precisely this term. G20 has K contributions as well.

The Hamiltonian constraint integrates as

    integral [R[h]+K²-K_ij K^ij]dvol_h=16 pi G E_m

on each exhaustion with convergent integrals (otherwise retain finite-domain
equality). There is no generally valid conversion of integral R[h] alone
into ADM mass with a useful sign. At K=0 it gives R[h]>=0; in a nonmaximal
slice it need not. The initial-data positive-mass theorem instead uses the
full constraints and dominant energy condition. It supplies M_ADM>=0 under
its global hypotheses, but G22 is compatible with that inequality. Thus
these coupled identities give necessary balances and no sign contradiction
for Omega!=0. They do not prove back-reaction supports localization.

## Final classification and decisive gap

| Result | Classification |
|---|---|
| General stationary scalar/stress/Einstein reduction, including nonzero K and shift | VERIFIED |
| Finite matter energy forces z in H¹ and L² for Omega!=0 | VERIFIED under G1 |
| Monopole coefficient -4G M_ADM Omega²/R in transformed V_eff | CONDITIONAL on G12/G14; algebra VERIFIED |
| Spherical exterior L² exclusion retaining Coulomb phase | VERIFIED under G2; CONDITIONAL geometric applicability |
| Global continuation and exclusion | OBSTRUCTION IDENTIFIED under G2's ellipticity and smoothness |
| Applying short-range Rellich directly to the gravitational monopole | OBSTRUCTION IDENTIFIED to that proposed proof |
| General anisotropic stationary coupled exclusion or existence | NOT YET ESTABLISHED |
| Universal stationary no-go beyond the rotational representation | NOT YET ESTABLISHED |
| Stability, independent review, novelty and PRD readiness | NOT YET ESTABLISHED |

Strongest defensible conclusion: gravitational back-reaction does not localize
a nonzero real-frequency rotational mode in G2's controlled class, despite
an attractive mass monopole. The original unrestricted stationary question
remains a precisely defined quadratic-pencil spectral and continuation
problem, coupled to (G2), (G7)–(G10). The general assumptions in the request
do not by themselves justify G12, spherical exterior coefficients or absence
of ergoregions. No existence inference is made from those gaps.

One decisive unresolved question: do the actual coupled smooth finite-energy
stationary Einstein–DEE solutions enforce an exterior coefficient class and
interior continuation property sufficient to exclude every positive-frequency
L² mode, without the extra spherical/strictly-timelike hypotheses of G2?

Conservative Physical Review D assessment: **not ready**. There is no general
stationary theorem, localized solution, stability analysis, independent
mathematical review or exhaustive priority assessment. This bounded extension
uses established ADM, Coulomb asymptotics and elliptic continuation methods.
No new general spectral theorem or independent peer review is claimed. No
Stage H, simulations, action modification or manuscript work is undertaken.

## Reproduction and observed verification record

Only the three authorized Stage G files are additions. Complete actual command
stdout/stderr and exit codes are appended below after execution; temporary logs
are outside the repository. A summary PASS is not an extra check. Historical
scope limits and the incomplete provenance-recovery gate remain unchanged.

### Observed final runs, 2026-11-10 UTC

Python 3.12.14, SymPy 1.14.0. All required commands exited 0; no FAIL occurred. The Stage G verifier had successful intermediate 59- and 61-check runs before the final 60-check version (one redundant algebra check was removed). The final count is 60 named checks, some covering multiple tensor components. Analytical global proofs are not counted.

Command: `python DEE_STAGE_G_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

```text
SymPy 1.14.0
PASS: general ADM inverse
PASS: ADM block congruence
PASS: unit triangular shear determinant
PASS: det g = -N squared det h via congruence
PASS: Lorentz kinetic contraction
PASS: coordinate Ttt
PASS: all coordinate Tti
PASS: positive ADM energy
PASS: all ADM momentum projections
PASS: all spatial stresses
PASS: active mass epsilon plus S
PASS: trace reversed spatial source
PASS: normal phase kinetic shift signs
PASS: Komar rotation charge cancellation
PASS: stationary stress common rotation
PASS: metric variation stress (0, 0)
PASS: metric variation stress (0, 1)
PASS: metric variation stress (0, 2)
PASS: metric variation stress (0, 3)
PASS: metric variation stress (1, 1)
PASS: metric variation stress (1, 2)
PASS: metric variation stress (1, 3)
PASS: metric variation stress (2, 2)
PASS: metric variation stress (2, 3)
PASS: metric variation stress (3, 3)
PASS: general 3D covariant phase density reduction
PASS: unexpanded phase current equivalence
PASS: amplitude action versus covariant reduction
PASS: flat Helmholtz limit
PASS: vacuum radial mass
PASS: amplitude perturbation polynomial
PASS: formal Hermitian Green boundary identity
PASS: amplitude boundary product rule
PASS: momentum integration product rule
PASS: phase real Green boundary and all shift terms
PASS: monopole volume density
PASS: monopole spatial density cancellation
PASS: monopole frequency density
PASS: isotropic optical frequency coefficient
PASS: effective Coulomb potential sign
PASS: exact radial density conjugation
PASS: spherical harmonic radial equation
PASS: Coulomb comparison residual 1
PASS: Coulomb comparison residual -1
PASS: Coulomb comparison Wronskian
PASS: asymptotic amplitude stiffness
PASS: nonmaximal K trace
PASS: Hamiltonian direct curvature Gauss
PASS: momentum direct Codazzi 0
PASS: momentum direct Codazzi 1
PASS: momentum direct Codazzi 2
PASS: stationary normal Ricci lapse with beta DK
PASS: spatial Ricci ADM (0, 0)
PASS: spatial Ricci ADM (0, 1)
PASS: spatial Ricci ADM (0, 2)
PASS: spatial Ricci ADM (1, 1)
PASS: spatial Ricci ADM (1, 2)
PASS: spatial Ricci ADM (2, 2)
PASS: lapse source coefficient
PASS: Hamiltonian matter coefficient
PASS: Stage G 60 exact local checks; global proofs are conditional as stated in the reports.
```

Command: `python DEE_STAGE_F_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

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

Command: `python DEE_STAGE_E_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

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

Command: `python DEE_STAGE_D_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

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

Command: `python DEE_STAGE_C_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

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

Command: `python DEE_STAGE_B_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

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

Command: `python DEE_STAGE_A_HALF_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

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

Command: `python DEE_STAGE_A_HALF_BOUNDED_VERIFIER.py`; exit code **0**. Complete stdout/stderr:

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

Command: `python evidence/einstein_dee_sources/verify_evidence.py --package-only`; exit code **0**. Complete stdout/stderr:

```text
Git preservation NOT CHECKED in package-only mode
PASS: 49 package hashes, 20 archive members, pinned source copies
Recovery gate: INCOMPLETE. Physical tests: NOT RUN; approval required.
```

The Stage C log contains 23 named exact checks; Stage B and the two Stage A½ verifiers report their original final PASS without invented numbered totals. Frozen provenance checks 49 package hashes, 20 archive members and pinned source copies. Its recovery gate remains INCOMPLETE; physical tests remain NOT RUN. Package-only mode explicitly does not perform Git preservation; the separate exhaustive comparison below does.

### Exact-parent content and mode preservation

Read-only comparison command: `python /tmp/stage_g_preserve.py`; exit code **0**. It enumerates `git ls-tree -rz` at Stage F, recomputes every worktree Git blob digest from its original bytes (symlink targets where applicable), checks modes, checks the index against the baseline, and compares the complete added-path set against the three authorized names. Actual output:

```json
{
  "baseline": "767783984a8502d892cf52048cd43eadde79e552",
  "baseline_tree": "0da0d813d76108bc880d2e34c47f82b261026cd1",
  "head": "767783984a8502d892cf52048cd43eadde79e552",
  "tracked_blobs_checked": 1120,
  "run43_named_blobs": 233,
  "content_or_mode_mismatches": [],
  "index_mismatches": [],
  "additions": [
    "DEE_STAGE_G_COUPLED_STATIONARY_EQUATIONS.md",
    "DEE_STAGE_G_SPECTRAL_LOCALIZATION.md",
    "DEE_STAGE_G_VERIFIER.py"
  ],
  "exactly_three_authorized_additions": true
}
```

All **1,120 Stage F tracked blobs** remain byte-for-byte unchanged, including **233 Run 43 named paths** and every other tracked artifact regardless of naming. Stage F historically compared 1,117 blobs against Stage E; the present baseline includes the three Stage F files and therefore has 1,120. This includes the original manuscripts/action, earlier reports/verifiers, archives, checkpoints, numerical evidence and historical outputs. No experiment was rerun.

Commands `git diff --check` and `git diff --exit-code 767783984a8502d892cf52048cd43eadde79e552 --` exited **0**, empty output before staging. Actual `git status --short` (exit **0**) before staging:

```text
?? DEE_STAGE_G_COUPLED_STATIONARY_EQUATIONS.md
?? DEE_STAGE_G_SPECTRAL_LOCALIZATION.md
?? DEE_STAGE_G_VERIFIER.py
```

Only these three Stage G paths are staged and committed. The exact-parent/index/content comparison is repeated before commit; remote commit and the three expected blob IDs are verified after push, and clean worktree status is checked. The final commit SHA and remote verification are reported in the delivery message because a commit cannot embed its own SHA. No merge or force push is authorized or performed. Stop after Stage G.
