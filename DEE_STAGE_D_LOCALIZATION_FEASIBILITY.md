# Stage D: localization feasibility and theorem limits

Date: 2026-11-10 UTC. Parent `514f5c77cb7250b064b54b3c989bde043c53ad70`.
Equations, conventions and axis/asymptotic conditions are in
`DEE_STAGE_D_RELAXED_SYMMETRY.md`. No action or prior evidence is changed.

**Scientific verdict:** Stage D proves a broader, explicitly conditional static
no-go and derives necessary conditions for other sectors. It produces no new
localized existence result. Relaxing translation invariance allows local current
redistribution, but does not evade the global static obstructions below.
A universal theorem covering the full Einstein–DEE theory is NOT YET ESTABLISHED.

## Theorem D1: phase rigidity and positive-domain amplitude rigidity

Let a connected static spatial domain have N>0 and positive definite M(rho).
Assume smooth symmetry-inheriting fields, globally single-valued real phase
lifts, and an exhaustion whose full boundary terms obey

    lim integral_boundary N (Q-Q_infty)^T M D_n Q dA = 0.

All boundaries are included; a regular filled axis has zero flux in cylindrical
coordinates. Circle-valued phases admit lifts on a simply connected slice;
nontrivial topology is an extra issue, not silently assumed away. Let the
phase kinetic integral exist. Integrate each phase equation
`D_i(N M D^i Q)=0`, multiply by `(Q-Q_infty)^T` and sum:

    integral N (D_i Q)^T M D^i Q dvol_h
      = lim integral_boundary N (Q-Q_infty)^T M D_n Q dA = 0.

The integrand is nonnegative and smooth; strict positivity gives DQ=0
everywhere, so Q is constant. Coupling through rho produces no missing volume
term: integration by parts acts on the complete current N M DQ. This proof
allows arbitrary smooth positive M, without a sign assumption on M'.
The verifier checks the exact local product-rule identity behind this integral.
Conservation alone without the boundary/lift assumptions does not imply rigidity.

Now take the original polynomial potential, a nonnegative amplitude rho>=0,
asymptotic rho=v, and also
`lim integral_boundary N(rho-v)D_n rho dA=0`, with convergent bulk integrals.
With DQ=0, multiply `D_i(N D^i rho)=N V'` by rho-v:

    integral N [|D rho|²+(rho-v)V'(rho)] dvol_h = 0,
    (rho-v)V'(rho)=lambda_V rho(rho+v)(rho-v)² >=0.

Thus rho is constant, and its asymptotic value makes it v. This proof does not
need Einstein's equations, convexity of V on the entire real line, or A>=0.
It applies on any prescribed static geometry satisfying its hypotheses. For a
signed scalar crossing negative values the displayed sign is not valid; no
such extension is made. Analogous reasoning for a negative-only field uses
the -v vacuum. The matter fields are vacuum, not a DEE transport equilibrium.

**VERIFIED / OBSTRUCTION IDENTIFIED:** under these boundary and field-domain
hypotheses, no nontrivial static localized matter solution exists even after
allowing both r and z dependence. This includes exact compact vacuum matching
on a globally regular isolated slice. Local bounded-domain boundary-driven
solutions are not excluded, because their boundary integral need not vanish.

## Theorem D2: gravitational static no-soliton obstruction

Assume a smooth complete orientable asymptotically Euclidean three-dimensional
slice with one end, no inner boundaries, the ADM hypotheses of the positive
mass theorem, a globally static Einstein solution with everywhere positive
lapse normalized to one at infinity, symmetry-inheriting scalar fields, M>0,
V>=0, and sufficient decay for the following integrals and the equality of
static Komar and ADM masses. These are substantive geometric assumptions.
They are not established by local matter continuation or by a symbolic program.

Write `C_ij=D_i rho D_j rho+D_i Q^T M D_j Q` and Kmat=trace_h C.
Then epsilon=Kmat/2+V and trace_h S=-Kmat/2-3V, hence

    epsilon+trace_h S=-2V,
    D²N=4pi G N(epsilon+trace_h S)=-8pi G N V.

Integrate the lapse equation on the complete slice. There are no horizon or
inner-boundary contributions. Its asymptotic surface flux gives

    M_ADM=M_K=(1/(4pi G)) lim integral D_n N dA
                =-2 integral N V dvol_h <=0.

The Hamiltonian constraint gives `R[h]=16pi G epsilon>=0`. The positive mass
theorem in dimension three therefore gives M_ADM>=0; its zero-mass rigidity
makes h Euclidean. The constraint then gives epsilon=0, so V=0 and all scalar
derivatives vanish by M>0. The normalized smooth lapse is one (also follows
from its harmonic equation and the asymptotic condition). Spacetime is
Minkowski with constant vacuum matter. QED.

This argument covers the signed polynomial scalar too, and does not require
global phase lifts: the stress and local derivatives suffice. It does not
assume axisymmetry or a special M' sign. It uses the established positive mass
theorem and static mass equality, rather than pretending to prove those global
geometric theorems symbolically. Static positive scalar kinetic energy and
V>=0 also obey the dominant energy condition; scalar curvature positivity is
the time-symmetric positive-mass input actually needed here.

**VERIFIED as a theorem with these hypotheses / OBSTRUCTION IDENTIFIED:** no
regular isolated static soliton exists in this sector. **CONDITIONAL:** its
application to any proposed geometry requires verifying every global hypothesis.
It excludes a larger class than Stage C and D1, but is a standard static
Einstein–scalar no-soliton argument specialized to the DEE sigma-model metric.
It is not a new theory-wide uniqueness theorem.

## Competing paths and what remains open

Path A has no nontrivial candidate within D1/D2. Flat-space vacuum perturbations
are not existence evidence for a DEE lump: near rho=v,
`(Delta-2lambda_V v²)delta rho=0`, `M(v)Delta delta Q=0` for static matter.
The same decaying regular energy identities eliminate these linearized modes.
A Yukawa tail `e^{-mR}/R` is an exterior possibility but singular at its center;
a tail alone supplies no interior. Time-dependent massive amplitude waves and
massless phase waves are ordinary scalar fluctuations, not stationary solitons.

The full stationary class is larger than the static metric selected here.
For example circular stationary coordinates can have
`W²(dphi-omega dt)²` instead of W²dphi². A nonconstant omega gives nonzero
extrinsic curvature and changes the Hamiltonian/momentum system; it was not
solved in Stage D. Matter with `Q=w t+q(r,z)` does not inherit time symmetry,
though its stress may. With polynomial M(v)>0 and asymptotically flat N->1,
a constant nonzero w leaves positive energy at infinity and cannot localize.
Likewise winding costs axis energy since M(0)=I. These two observations exclude
simple borrowed boson-star/winding recipes, not every possible stationary
symmetry realization. No claim that all stationary metrics are static is made.

Path B therefore succeeds in the isolated static sector, with two explicit
proofs. Horizons, physical boundaries, non-asymptotically-flat ends, singular
maps, non-inheriting fields, topology in D1, negative potentials or nonpositive
kinetic sectors change hypotheses. They are neither proposed physical solutions
nor modifications authorized to rescue localization. D2 still covers topology
when its positive-mass and completeness hypotheses hold. Unknown original
coupling functions with negative V cannot inherit the V>=0 conclusion by fiat.

For initial data alone, localized smooth matter with zero momenta has j=0,
but its Hamiltonian constraint still requires an actual elliptic solution.
Constraint compatibility is not equilibrium: D2 precludes its static completion
in the stated class. Eventual second-order or nonlinear verification would need
all Einstein constraints, lapse/spatial equations, both phase EOM, amplitude
EOM, controlled residuals and boundary limits at the claimed order, positive
kinetic eigenvalues, regular geometry and, separately, a stability analysis.
The Stage A½ and B leading response solves none of those global problems.
No nonlinear evolution or large parameter/PDE search is performed here.

## Relation to established results and publication significance

The field-space metric is `d rho²+dQ^T M(rho)dQ`: this is an Einstein nonlinear
sigma model with potential. Off-shell conservation, phase Green identities and
the static active-mass cancellation are standard structures, not new DEE laws.
In flat three-space scaling `fields_lambda(x)=fields(x/lambda)` gives
`E(lambda)=lambda E_gradient+lambda³ E_potential`; stationarity would require
`E_gradient+3E_potential=0`. For nonnegative terms this is Derrick's obstruction,
assuming an admissible differentiable scaling variation and finite energy.
It is not used as a curved-spacetime proof; D2 supplies that argument.

Relevant comparison anchors:

* G. H. Derrick, J. Math. Phys. 5, 1252 (1964),
  https://doi.org/10.1063/1.1704233 — flat-space scaling obstruction.
* D. J. Kaup, Phys. Rev. 172, 1331 (1968),
  https://doi.org/10.1103/PhysRev.172.1331; R. Ruffini and S. Bonazzola,
  Phys. Rev. 187, 1767 (1969), https://doi.org/10.1103/PhysRev.187.1767 —
  scalar gravitational configurations with time-dependent field structure.
* R. Schoen and S.-T. Yau, Commun. Math. Phys. 65, 45–76 (1979),
  https://doi.org/10.1007/BF01940959 — positive mass input;
  E. Witten, Commun. Math. Phys. 80, 381–402 (1981),
  https://doi.org/10.1007/BF01208277 — alternative positive-energy proof.

Boson stars use complex scalars whose time-dependent phase can have vanishing
stiffness as their amplitude tends to zero. Here M(0)=I and M(v)>0, so a
spatially constant temporal frequency persists energetically at infinity.
This concrete difference, vacuum at v rather than zero, and extension of the
current integral to the mixed M are model-specific applications, not established
novel mathematics. No exhaustive priority search or independent peer review
has been performed; originality is NOT YET ESTABLISHED. No unsupported named
uniqueness theorem is applied beyond its hypotheses.

A narrowly focused PRD manuscript is not ready on these results alone: the
negative static result is scientifically useful but largely a specialization
of established arguments. A credible research route would require a careful
literature/priority audit and a sharply defined unresolved sector of the
unchanged action, with an actual existence or exclusion theorem. Stage D does
not authorize submission or start that work.

## Verification, classification and stopping point

`DEE_STAGE_D_VERIFIER.py` uses exact SymPy algebra. It derives all three reduced
Euler–Lagrange equations independently from the action density, checks all four
off-shell covariant stress divergences using Christoffel symbols, both phase
currents, mixed scalar derivatives, mixed Einstein geometry, Hamiltonian and
momentum geometry, lapse curvature, integral product rules, potential signs as
factorizations, vacuum values, determinant and perturbative mass. Inequalities,
global lifts, boundary limits, ADM/Komar equality, completeness and positive-mass
rigidity are analytical hypotheses/theorem inputs, not machine-certified facts.

| Conclusion | Classification |
|---|---|
| Full relaxed static equations and exact identity comparisons | VERIFIED |
| Local axial flux redistribution allowed by conservation and closure | VERIFIED; no localized solution follows |
| Phase rigidity with vanishing full boundary term and lifts | VERIFIED / OBSTRUCTION IDENTIFIED under D1 |
| Nonnegative-amplitude polynomial static rigidity | VERIFIED / OBSTRUCTION IDENTIFIED under D1 |
| Complete isolated static Einstein no-soliton result | VERIFIED / OBSTRUCTION IDENTIFIED under D2 |
| A proposed metric satisfies D1/D2 global assumptions | CONDITIONAL; must be checked for that proposal |
| Nontrivial localized matter, solved Einstein initial data, full spacetime | NOT YET ESTABLISHED |
| Universal stationary/dynamical no-go, stability, originality | NOT YET ESTABLISHED |

Recommended next scientific test, only after separate authorization: analyze
symmetry inheritance for stationary stress in this positive-kinetic model and
whether the permitted asymptotic vacuum forces temporal phase gradients and
rotation to vanish. This is a small analytical test that would determine whether
D2 can be extended, before considering any PDE construction. Stop after Stage D.

Observed verification, 2026-11-10 UTC (Python 3, SymPy 1.14.0):

* Stage D: exit 0, all 27 exact checks pass.
* Unchanged Stage C, Stage B, Stage A½ symbolic and Stage A½ bounded verifiers:
  each exit 0, final PASS. These are regressions, not nonlinear verification.
* Source provenance verifier with `--package-only`: exit 0, 49 package hashes,
  20 archive members and pinned source copies pass. Its full historical mode
  exits 1 at `addition outside evidence package`: that guard forbids later
  research additions already present at the Stage C parent. This is a scope
  incompatibility, not a hash failure; the verifier is preserved unchanged.
* `git diff --check`: pass. Exact diff against the Stage C parent, excluding only
  the three named new Stage D files: empty. This separately verifies preservation
  of all baseline files, including Run 43 and frozen evidence.

Reproduce with `python DEE_STAGE_D_VERIFIER.py`, the four prior verifiers named
above, and `python evidence/einstein_dee_sources/verify_evidence.py --package-only`.
The historical recovery gate printed by that last program remains unchanged;
no missing historical provenance or physical test is claimed recovered here.
