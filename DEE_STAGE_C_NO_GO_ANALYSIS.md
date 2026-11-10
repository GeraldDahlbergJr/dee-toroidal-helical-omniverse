# Stage C: no-go theorems and limits

Date: 2026-11-10 UTC. Base: `0a9c4f84c383089d3f74b49fe2a4d5fb05c75856`.
The action, notation, derivation and source authority are in
`DEE_STAGE_C_VACUUM_TRANSITION.md`. These are analytical theorems with explicit
hypotheses; symbolic verification checks identities, not global existence.

## Theorem 1: vacuum rigidity in the positive kinetic sector

Let a C² solution of the original action on a connected Lorentzian open set
have positive definite M(rho) everywhere. If T_mu nu=0 there, then rho and
both phases are locally constant, V=V'=0. For circle-valued phases, constant
means a constant map, independent of the choice of local phase representative.
Conversely these conditions give a matter vacuum.

**Proof.** For every null l,
T(l,l)=(l rho)²+(l Q)^T M(l Q). Take l=e_0+e_i and e_0-e_i
for i=1,2,3 in a local orthonormal frame. Vanishing and positivity imply
d_0 rho +/- d_i rho=0 and d_0 Q +/- d_i Q=0. Thus all derivatives
vanish. Connectedness makes the fields constant. T=-Vg gives V=0 and
the amplitude EOM gives V'=0. The converse follows by substitution. This
proof does not assume a static geometry, V>=0, or a diagonal metric. QED.

## Theorem 2: finite axial transition obstruction

Assume a connected smooth region contains a radial interval reaching a finite
interface and an axial translation generator xi=partial_z, nonvanishing in
regular coordinates there. Assume dQ is invariant under xi and all transverse
coordinate translations used in the radial ansatz, and Q has the retained
form w t+m phi+k z+q(r). Assume k!=0 and M positive definite through the
interface and prospective open vacuum exterior. No continuous-gradient,
source-free finite transition in this family reaches vacuum.

**Proof.** The constants w,m,k cannot be radial functions: for example
(d dQ)_rz=partial_r(d_z Q)-partial_z(d_r Q)=k'(r)=0 by invariant
gradients. Equivalently Cartan's identity gives d(xi Q)=0 when Lie_xi dQ=0.
Hence dQ(xi)=k stays nonzero throughout. Vacuum rigidity requires dQ=0
on the exterior and, by continuous traces, at its finite interface. This
contradicts k!=0. The argument is independent of V, amplitude response,
Einstein gauge and collar width. In a finite regular diagonal metric a
positive lower energy bound, derived in the companion file, supplies a second
proof. QED.

This theorem applies directly to the established non-null axial matter family
and any radial continuation retaining its gradients. It also excludes the
common null axial family as a vacuum, even though its kinetic contraction
vanishes. It does not require extrapolating the Stage B weak-field metric.
A discontinuous phase cutoff is not an exception: it introduces singular
field derivatives and lies outside the hypotheses.

## Theorem 3: no smooth compact radial matter support in the static sector

Assume the diagonal static radial metric of the companion file is smooth and
nondegenerate near a finite R>0, with positive N,H,L,Z. Let V,M be C²,
M invertible near rho_v, and V(rho_v)=V'(rho_v)=0. Consider a classical
radial solution with invariant phase gradients that matches smoothly to
constant vacuum fields for r>=R, without scalar surface sources. Then that
solution is vacuum on its connected regular radial interval (continued only
while these hypotheses hold).

**Proof.** The prescribed constant exterior phases and phase-trace continuity
set w=m=k=0. Normal
phase-current continuity fixes j=0 and invertibility gives q'=0. At R the
amplitude has rho=rho_v and rho'=0. Its equation is the first-order system

    rho'=p,
    p'=-(log(S/H²))' p + H² V'(rho).

The coefficients are regular and the right-hand side locally Lipschitz in
(rho,p). The constant solution (rho_v,0) has the same Cauchy data; ODE
uniqueness gives equality on a neighborhood, then on the connected regular
interval by continuation. Treating the metric of any hypothetical coupled
solution as given smooth coefficient functions is legitimate: the constant
amplitude solves this equation for that same metric. Its T is zero, so any
coupled solution there must also have a vacuum metric. QED.

This is stronger than a failure to find a radial shooting solution, but is
limited to a regular static radial system and exact finite vacuum matching.
It does not exclude tails that only approach vacuum at infinity, dynamical
compact-support initial data, nonradial fields, characteristic boundaries,
non-Lipschitz potentials, degenerating metrics, or shells. It is not a
unique-continuation theorem for arbitrary Einstein–matter PDEs.

## Currents, energy and remaining mathematical questions

For a source-free radial solution, j=S M q'/H² is constant. A vacuum
exterior sets it to zero; a nonzero flux cannot simply disappear in a collar.
In general, integrate nabla_mu J^mu=0 over a spacetime volume: incoming
current must be balanced by outgoing current or changing charge. This does
not require all axial currents of a different z-dependent solution to stay
constant pointwise. Radial or time-dependent redistribution is possible in
principle and has not been solved. Shift charge conservation alone is not a
universal localization no-go.

For static z-independent matter in the positive-energy sector, infinite z
extent precludes finite total matter energy unless the energy per length
vanishes. For nonzero winding on the full disk, phase topology and M>=I
independently obstruct regular finite axis energy. Neither argument is a
proof excluding every nonwinding, z-dependent, time-dependent solution.

A decisive result beyond these theorems needs a specified symmetry class for
fields versus stress, an admissible field domain and topology, fixed coupling
functions/parameters, actual constraint-satisfying Einstein interior data,
outer geometry and falloff, and boundary/Cauchy data for all phases and rho.
To prove existence requires a boundary-value/PDE existence argument with
kinetic, geometric and regularity control throughout the transition. To prove
a broader static no-go would require applicable global hypotheses (for example
completeness, lapse behavior, no inner boundaries and asymptotic conditions)
and an actual integral or uniqueness proof. Those conditions are not supplied
by local ODE continuation or the Run 43 evidence. No such theorem is invented.

| Finding | Status and scope |
|---|---|
| Exact vacuum conditions | Proved for strictly positive M |
| Amplitude taper with retained nonzero k | Impossible at finite regular interface |
| Smooth compact radial static matter support | Impossible under Theorem 3 hypotheses |
| T approaching zero at infinite radius | Excluded with uniform kinetic bound and bounded Z for retained k; otherwise unresolved |
| Phase taper k(r)z in current ansatz | Violates invariant gradients and generates z-dependent energy |
| General Einstein–DEE transition after symmetry change | Unresolved; no theory-wide obstruction established |
| New physical exterior or complete spacetime | Not constructed |

## Reproducible verification and evidence

Run from the repository root with Python 3 and SymPy 1.14.0:

    python DEE_STAGE_C_VERIFIER.py
    python DEE_STAGE_B_VERIFIER.py
    git diff --check

The Stage C verifier uses exact symbolic expressions only. It checks reduced
Euler–Lagrange equations against the covariantly derived reduction, all four
off-shell stress divergences by direct connection calculations, phase-gradient
integrability, the vacuum/axial polynomial identities, the ODE vacuum fixed
point, and the normal Einstein constraint by independent Ricci and extrinsic-
curvature calculations. Positivity, topology and uniqueness are proved in
text, not inferred from symbolic samples. The prior Stage B verifier is a
compatibility regression only and is not evidence for an exact Einstein
interior. No numerical PDE, full evolution or Run 43 execution is performed.
Only the three Stage C deliverables are added; all earlier files are preserved.
Stage C stops here. No PRD readiness or next-stage authorization is claimed.

Observed verification (2026-11-10 UTC): Stage C exited 0 with SymPy 1.14.0
and 23 exact identity checks passing; the unchanged Stage B verifier exited 0
with all 26 checks passing. `git diff --check` passed. An initial verifier
implementation expanded an odd polynomial before substituting its branch
relation, preventing the intended substitution; the expression was corrected
to substitute in the even bracket. No physical equation was changed.
