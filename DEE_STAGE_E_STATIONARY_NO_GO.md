# Stage E: necessary conditions, rotation and localization obstruction

Date: 2026-11-10 UTC. Exact parent:
`570a9fdcf8ce7fb62e8d377bb535e6da6b42595d`.
Branch: `research/stage-e-stationary-analysis`.
Conventions and full equations: `DEE_STAGE_E_STATIONARY_ANALYSIS.md`.

**Scientific result — OBSTRUCTION IDENTIFIED, with explicit hypotheses:** the
regular, filled-axis, horizonless, asymptotically flat circular stationary class
with rho=rho(x^A), Q=w t+m phi+q(x^A), strictly positive M and the global
conditions below contains only Minkowski spacetime with constant vacuum matter.
The temporal and angular obstructions are independent of the static no-soliton
argument. Staticity follows from a separate rotation integral, not from zero
matter angular momentum alone. This is a conditional extension of Stage D,
not a universal stationary Einstein–DEE no-go theorem.

Throughout, **VERIFIED** denotes an analytical proof under stated hypotheses
and/or a specified exact local identity check; it does not mean external peer
review. **CONDITIONAL** marks application-dependent assumptions.
**OBSTRUCTION IDENTIFIED** marks exclusion within those assumptions.
**NOT YET ESTABLISHED** marks a question left open, never positive evidence of
existence. The unchanged polynomial action, signature and field definitions
are specified in the analysis report. No earlier research or Run 43 is changed.

## Theorem E1: temporal-frequency necessary condition

**Statement (VERIFIED / OBSTRUCTION IDENTIFIED).** Let the stationary field
representation and circular metric of the analysis report hold on an
asymptotically Euclidean end, with constant finite w,m, N->1, omega->0,
rho->v, and M(rho)->M(v)>0. Assume the original V>=0 and h positive. If the
observer energy epsilon=T(n,n) tends to zero at infinity, then w=0. More
strongly, if w!=0, epsilon has a strictly positive lower bound outside a
compact set, and the proper matter energy on that end is infinite. Thus
asymptotically flat vacuum stress and finite matter energy in this class
require w=0.

**Proof.** In an orthonormal frame based on the slice normal,

    epsilon = V + 1/2 [gamma^AB rho_A rho_B
                  +gamma^AB q_A^T M q_B
                  +m^T M m/W²+U^T M U/N²], U=w+omega m.

All terms are nonnegative. The last quadratic converges to w^T M(v)w.
If w!=0 this limit is strictly positive, so, for example,
`epsilon >= w^T M(v)w/4` sufficiently far out. The end has infinite volume.
This proves both conclusions without assuming the other gradients vanish.
If spatial gradients and V also vanish, the full orthonormal limiting tensor
is explicitly

    T_hat a hat b -> diag(B/2,B/2,B/2,B/2), B=w^T M(v)w>0,

with zero limiting mixed entries when their gradients decay. Thus this is
nonzero stiff scalar stress, not a cosmological constant or vacuum. Standard
differentiated asymptotic Einstein falloff also excludes a nonzero limiting
stress. Finite ADM mass alone, without the matter/metric falloff needed for
these conclusions, is not used as an unexplained integral-energy identity.
QED.

A zero Lorentz kinetic contraction is irrelevant: in flat spacetime take
rho constant and d_t Q=d_z Q=k. Then the phase kinetic scalar is zero but
T_tt=T_tz=T_zz=k^T M k (with V=0). Positivity and the full tensor prohibit
calling such a null gradient vacuum. The verifier checks these components.
A constant angular velocity at infinity could change U's coordinate expression,
but the asymptotic time normalization here fixes omega->0. It is not an escape
from this theorem in its specified inertial frame.

**Boundary of E1 (CONDITIONAL / NOT YET ESTABLISHED).** If the time dependence
is not an additive constant phase frequency, E1 alone does not apply. Neither
degenerate M(v), nonvacuum ends nor other asymptotic normalizations have been
proved impossible here; they change hypotheses. No kinetic function is modified.

## Theorem E2: angular winding at a regular filled axis

**Statement (VERIFIED / OBSTRUCTION IDENTIFIED).** Suppose a smooth rotation
axis is filled, the spacetime and rho extend regularly to it, N is positive
and finite there, M(rho) is continuous and strictly positive there, and
Q=w t+m phi+q(x^A) defines globally regular real scalars or regular maps into
a product of circles. Then m=0. Independently, m!=0 has divergent local proper
matter energy near any nonzero-length regular axis segment, even if one only
requires the gradient stress on the punctured neighborhood.

**Proof for real scalars.** At fixed t,x^A, single-valuedness under
phi->phi+2pi requires 2pi m=0 component by component. Real scalars cannot
have a nonzero constant winding even on an annulus without permitting a
multivalued scalar (which is a different field definition).

**Proof for circle maps.** In the usual 2pi phase normalization m is an
integer two-vector. An axial orbit bounds a disk through the filled axis.
A continuous circle map on that disk restricts to degree zero on its boundary;
therefore both degrees vanish. Other phase periods merely rescale this statement.
Alternatively, the axial Killing field vanishes at the axis, so a smooth scalar
map has dQ(eta)=0 there; the ansatz keeps this derivative equal to m.

**Independent energy proof.** Around a regular axis segment, proper transverse
distance s obeys W=s+O(s³). Smooth rho and strict positivity give a uniform
local lower eigenvalue bound M>=mu I with mu>0. Hence

    epsilon >= m^T M m/(2W²) >= c |m|²/s².

The proper volume is proportional to s ds dphi dl, so a nonzero m gives
`integral_0^delta ds/s=infinity`. Neither temporal terms nor rotation cancels
this positive contribution. Letting rho vanish on the axis does not help:
M(0)=I. QED.

**Additional asymptotic obstruction (VERIFIED under the stated end bounds).**
Even if the axis is excised, a constant nonzero m on a full asymptotically
Euclidean end has divergent total matter energy whenever the end contains
an angular belt of fixed positive solid angle. On that belt W is comparable
to R, M>=mu I, epsilon>=c |m|²/R² and dvol_h is comparable to R² dR dOmega.
The radial integral diverges linearly. Stress may tend to zero there while
its total energy diverges. This is a distinct asymptotic obstruction, not
an argument that every circle map on every nontrivial topology is forbidden.

**Limits (CONDITIONAL / NOT YET ESTABLISHED).** An excised-axis annulus supports
smooth circle-valued winding locally; on a bounded annulus its energy can be
finite. This is a kinematic boundary example, not a self-gravitating isolated
DEE solution. If loops are noncontractible, the degree proof fails. Compact
cycles in a different interior topology need not behave as m dphi on a full
Euclidean end. Boundaries, defects and different ends must be analyzed with
their own current, Einstein and energy conditions. Strict positivity only at
each punctured point, with eigenvalues tending to zero at a missing axis,
would not supply the local uniform bound used here. That is not the smooth
filled-axis sector with the unchanged M and regular rho.

## Theorem E3: source-free circular rotation rigidity

**Statement (VERIFIED / OBSTRUCTION IDENTIFIED).** Let a connected smooth
circular stationary solution have a global spacelike slice with N>0,
h=gamma+W²dphi², smooth filled axes, and a global bounded smooth omega near
them. Suppose j_phi=0, no physical inner boundaries, and an exhaustion of the
meridional orbit space for which

    lim integral_boundary omega (W³/N) D_n^gamma omega dl_gamma = 0.

Assume the weighted gradient integral exists (or use the nonnegative exhaustion
identity with its boundary limit). Then D omega=0, and omega->0 on the end
implies omega=0 throughout. In particular w=m=0 implies j_phi=0 and satisfies
the source hypothesis; the isolated falloffs in the analysis report imply
the boundary condition.

**Proof.** The phi momentum equation is

    partial_A(sigma W³ gamma^AB omega_B/N)=-16 pi G sigma W j_phi.

Multiply by omega and integrate by parts on the full quotient:

    integral (W³/N)|D_gamma omega|² dvol_gamma
      = lim integral_boundary omega (W³/N) D_n^gamma omega dl_gamma
                     +16 pi G integral W omega j_phi dvol_gamma.

The right side vanishes under the hypotheses. N>0 and W>0 off the axis make
the left integrand nonnegative and force D omega=0. Smooth extension and
connectedness give a single constant, fixed to zero by the end condition.
There is no physical boundary contribution from a smooth filled axis: removing
a tube s=delta for the calculation gives W³=O(delta³), bounded N^-1 and
omega, and bounded derivatives (in fact transverse omega_s=O(delta)). The
term tends to zero. At infinity, W=O(R), omega=O(R^-3), Domega=O(R^-4),
dl_gamma=O(R dtheta), so the integrand integrated there is O(R^-3).
The weighted bulk is integrable at both ends. Include every other boundary
component in the exhaustion; the theorem excludes any unaccounted one. QED.

An invariant check uses beta=-omega eta and
`K_ij=(D_i beta_j+D_j beta_i)/(2N)`. Multiply the momentum constraint by beta^i:

    integral N K_ij K^ij dvol_h
       = integral_boundary beta^i K^j_i s_j dA
                               -8 pi G integral beta^i j_i dvol_h.

For purely axial beta the source is -omega j_phi; the same smoothness and
falloff kill the boundary. This reproduces E3 without assuming a single global
meridional chart. On quotient patches the fluxes on artificial interfaces
cancel. It requires a global circular shift and slice, as stated.

**Useful stronger identity (VERIFIED, conditional boundaries).** If w=0 but
m is retained, U=omega m and the rotation equation gives

    integral [(W³/N)|D_gamma omega|²
             +16 pi G (W/N) omega² m^T M m] dvol_gamma
        = lim integral_boundary omega (W³/N) D_n^gamma omega dl_gamma.

Thus its zero-boundary conclusion remains coercive with positive M. This
does not cure the filled-axis winding singularity, or remove boundary terms
on a domain where the axis is excised. No such domain is silently filled.

## Four distinct angular statements

**VERIFIED distinctions:** the following implications require their stated
extra inputs; no pair is interchangeable by definition.

| Statement | What it establishes and what else is needed |
|---|---|
| omega->0 at infinity | Fixes asymptotic rotation gauge; does not set Domega=0 |
| j_phi=0 | Makes the rotation equation homogeneous; global boundary data remain |
| J_ADM=0 | Fixes one unweighted asymptotic flux; does not alone kill the weighted bulk integral |
| omega=0 everywhere | Together with the chosen circular form and N>0 gives a static metric |

For w=m=0, the matter formula proves j_i=0 for all i. Integrating the angular
momentum constraint on a complete slice without other boundary fluxes then
gives J_ADM=0. The *separate* weighted E3 identity proves omega=0. One may
also derive J_ADM=0 after staticity; it is not used as a sufficient premise.
A merely vanishing unweighted surface integral of Domega would not force
`integral omega Domega` to vanish on an arbitrary inner boundary.

**Counterexample boundary (VERIFIED):** Kerr vacuum has j_phi=0 and omega->0
at infinity, but nonzero rotation and angular momentum. It solves the original
DEE theory with constant rho=v,Q as vacuum matter. Its exterior has a horizon
inner boundary and its full spacetime is not a regular horizonless soliton.
The horizon flux in the rotation identity cannot be discarded; N also tends
to zero there. Kerr refutes an unrestricted inference from zero matter source
and asymptotic nonrotation to staticity. It is not a counterexample to E3.
A rigidly rotating Minkowski coordinate chart has constant nonzero omega and
Domega=0, but violates the chosen asymptotic inertial normalization. These
examples specify exactly which hypotheses they change.

## Theorem E4: conditional stationary no-soliton extension of Stage D

**Statement (VERIFIED / OBSTRUCTION IDENTIFIED).** Consider a smooth solution
of the unchanged polynomial Einstein–DEE action with M(rho)>0 everywhere,
lambda_V,v,G>0, with the selected circular stationary metric and matter ansatz.
Assume:

1. A connected complete smooth orientable spacelike slice with one asymptotically
   Euclidean end satisfying the three-dimensional positive mass theorem
   hypotheses; no singularities, horizons or physical inner boundaries.
2. A global circular stationary decomposition with N>0, h positive and smooth,
   filled smooth elementary-flat axes, globally regular scalar maps (real or
   circle valued), and globally defined constant w,m as in the ansatz.
3. rho->v, M->M(v)>0, N->1, omega->0, asymptotic vacuum-frame stress and finite
   proper matter energy; the differentiated metric/scalar falloffs and integral
   bounds specified in the analysis report. In particular the rotation boundary
   term of E3 vanishes, and after staticity Komar and ADM masses agree with
   lapse flux `(4 pi G)^-1 lim integral D_s N dA`.

Then w=m=0, omega=0, the matter is constant vacuum and the spacetime is
Minkowski (up to coordinates). Consequently no nontrivial finite-energy
stationary soliton exists within this class. The proof does not require rho>=0,
A>=0, a sign of M', global real lifts of circle phases or simply connected
phase topology.

**Proof.** E1 sets w=0 and E2 sets m=0. Hence j_phi=0. E3 sets omega=0,
so beta=K_ij=0 and n rho=n Q=0. Only now are the Stage D static equations
available. Write X=h^ij C_ij>=0. They give

    epsilon=X/2+V, epsilon+S=-2V,
    R[h]=16 pi G epsilon>=0,
    D²N=-8 pi G N V.

There are no inner fluxes. Integrate the lapse equation, using the asymptotic
normalization and mass equality:

    M_ADM=M_K=-2 integral N V dvol_h <=0.

The three-dimensional positive mass theorem gives M_ADM>=0; its zero-mass
rigidity gives h Euclidean. The Hamiltonian constraint now gives epsilon=0,
so X=V=0, and strict M positivity sets every scalar derivative to zero. The
asymptotic value fixes rho=v; phases are constant maps. N is smooth harmonic
on Euclidean space and tends to one; the maximum principle gives N=1.
Therefore spacetime is Minkowski. QED.

Positive mass rigidity is an established global geometric theorem, not a
symbolic identity and not newly proved here. The stated complete three-dimensional
orientable smooth asymptotically Euclidean manifold with R>=0 is its standard
setting (with the usual decay/integrability conditions); no unverified
application to an explicit candidate geometry is made. The mass equality is
included with the differentiated stationary/static asymptotic hypotheses,
as in Stage D. Matter positivity ensures R>=0 only after the actual static
reduction in this proof, though the general constraint also gives
R[h]=K2+16 pi G epsilon>=0 for this maximal slice.

**Exact obstruction to a naive broader extension (VERIFIED).** Before E3,
the correct stationary lapse equation integrates to

    (4 pi G)^-1 lim integral D_s N dA
      =(4 pi G)^-1 integral N K2 dvol_h
                         +2 integral N(P-V) dvol_h.

Both K2 and P contribute positively and can defeat the static nonpositive
mass estimate. In a general stationary Komar decomposition the shift/momentum
contribution must also be retained. A noncircular shift can have other K
components, and stationary stress does not by itself imply nQ=0. Thus reusing
the Stage D static lapse equation would omit precisely the terms that require
proof here. In E4 they vanish by E1–E3; outside E4 they are an identified
mathematical obstruction to that proof, not a demonstration of existence.

## What remains unresolved

**NOT YET ESTABLISHED:** stationary stress in the full nonlinear sigma model
need not be assumed to arise from invariant phase gradients and invariant rho.
A field-space isometry can generate time dependence with invariant stress;
nonadditive symmetry realizations or a failure of circularity are not classified
by the selected ansatz. The kinetic metric is
`d rho²+dQ^T(I+rho² A)dQ`. Its always-present Q translations yield the w
considered here, but special couplings can give other isometries. At A=0,
for example, the phase plane is flat and has rotations: time-periodic Q can
have stationary stress without Q=w t+q. This is a symmetry example, not a
regular localized solution of the full equations. For real noncompact Q,
rotate its Cartesian components; for a torus phase target, such continuous
rotations need not be global isometries. The target and field topology matter.
Neither asymptotic vacuum nor positivity alone has been proved here to force
all such matter to inherit the stationary symmetry.

Noncircular stationary geometries, absence of a global positive-lapse slicing,
horizons, physical boundaries and other ends require separate arguments. Some
already have vacuum solutions (Kerr); none supplies an isolated regular matter
soliton by this observation. Nontrivial topology of q alone does not evade E4
when all its geometric and regularity hypotheses hold, since its gravitational
proof uses no phase lifts. Changing to negative V or an indefinite/degenerate
kinetic sector is outside this authorized positive sector, not a proposed fix.
No nonlinear solution, stability theorem or universal impossibility result is
established. A decisive E4 theorem has been obtained for the chosen class;
beyond that class the present analysis stops at these explicit obstructions.

**Single most important unresolved scientific question:** does every globally
regular, horizonless, asymptotically flat finite-energy stationary solution of
the original positive-kinetic theory admit the symmetry inheritance and global
circular reduction required by E4, or can a different stationary-stress scalar
realization evade them? This question must remain open here. No Stage F is begun.

## Relation to established results and publication significance

**VERIFIED comparison:** this is an Einstein nonlinear sigma model with a
nonnegative potential. The stress formula, ADM equations, Green identities,
axis degree argument and positive-mass input are standard. The circularity
conditions are not a theorem that all matter models have circular metrics:
orthogonal transitivity requires Ricci/matter compatibility and global
conditions. The mixed equations in this work explicitly expose those restrictions.
No vacuum/electrovacuum circularity or black-hole uniqueness result is applied
to DEE without its matter hypotheses.

Lichnerowicz-type stationary no-soliton arguments normally assume a globally
timelike stationary Killing field and suitable global and matter conditions.
A scalar potential can violate the strong energy condition, and stationary
sigma-model fields need not inherit the symmetry. Consequently a named
stationary-vacuum theorem is not a substitute for E1–E4. Here the local metric
may initially allow ergoregions; the source-free circular rotation identity
removes rotation under its own explicit hypotheses. Stage D's mass argument
is then used in its valid static setting. This supplies a direct proof rather
than claiming a stronger established theorem has automatically covered DEE.

Boson stars provide the relevant warning about stationary stress versus
stationary fields. A complex scalar written as f exp(i(Omega t+n phi)) has
phase stiffness proportional to f². As f->0 at spatial infinity, its frequency
can leave zero asymptotic stress; at an axis f can vanish sufficiently rapidly
to regularize a winding map of the full complex scalar. The polar phase is
then not itself a regular independent circle field at f=0. DEE instead has
two independent phase kinetic coefficients with M(0)=I and M(v)>0. It cannot
hide the additive w or filled-axis m using an amplitude zero. This distinction
does not rule out all other time dependence of the original three fields.
Stationary boson stars do not contradict E4's field representation and
asymptotic hypotheses.

Comparison sources (not an exhaustive priority search):

* A. Lichnerowicz, *Théories relativistes de la gravitation et de
  l'électromagnétisme*, Masson (1955): classical stationary no-soliton context.
* R. M. Wald, *General Relativity*, University of Chicago Press (1984),
  stationary/axisymmetric geometry, ADM/Komar charges and circularity framework.
* P. T. Chruściel, J. L. Costa and M. Heusler, *Stationary Black Holes:
  Uniqueness and Beyond*, Living Rev. Relativity 15, 7 (2012),
  https://doi.org/10.12942/lrr-2012-7,
  https://arxiv.org/abs/1205.6112: circularity and stationary uniqueness
  hypotheses, not an authorization to import electrovacuum results to DEE.
* D. J. Kaup, Phys. Rev. 172, 1331 (1968),
  https://doi.org/10.1103/PhysRev.172.1331; R. Ruffini and S. Bonazzola,
  Phys. Rev. 187, 1767 (1969), https://doi.org/10.1103/PhysRev.187.1767.
* F. E. Schunck and E. W. Mielke, *General relativistic boson stars*,
  Class. Quantum Grav. 20, R301–R356 (2003),
  https://doi.org/10.1088/0264-9381/20/20/201,
  https://arxiv.org/abs/0801.0307: stationary stress and rotating boson stars.
* R. Schoen and S.-T. Yau, Commun. Math. Phys. 65, 45–76 (1979),
  https://doi.org/10.1007/BF01940959; E. Witten, Commun. Math. Phys. 80,
  381–402 (1981), https://doi.org/10.1007/BF01208277: positive mass and rigidity.

Bibliographic metadata for the two linked arXiv reviews was checked in this
session; this is not a claim of an exhaustive contemporary literature review.
The proof itself is self-contained except for its explicitly identified
positive-mass and stationary/static mass-equality inputs.

**NOT YET ESTABLISHED (originality/publication):** applying these standard
structures to the mixed DEE matrix, identifying the unavoidable vacuum phase
stiffness, and chaining frequency/axis/rotation obstructions are model-specific
results. No new general mathematical theorem or priority is claimed. Symbolic
PASS results certify local equalities, not physical novelty, independent review,
stability, observational relevance or Physical Review D readiness. No PRD
manuscript is drafted or submitted.

## Reproduction, findings and preservation

From the repository root, Python 3 and SymPy 1.14.0:

    python DEE_STAGE_E_VERIFIER.py
    python DEE_STAGE_D_VERIFIER.py
    python DEE_STAGE_C_VERIFIER.py
    python DEE_STAGE_B_VERIFIER.py
    python DEE_STAGE_A_HALF_VERIFIER.py
    python DEE_STAGE_A_HALF_BOUNDED_VERIFIER.py
    python evidence/einstein_dee_sources/verify_evidence.py --package-only
    git diff --check

The unchanged provenance verifier's full historical mode has a known later-file
scope guard documented in Stage D; use its package-only mode for frozen hashes.
The direct baseline tree comparison separately checks every earlier tracked
file, including Run 43. No Run 43 experiment is rerun. Logs are execution
artifacts outside the repository; only the three Stage E deliverables are added.

| Finding | Classification |
|---|---|
| Stationary equations, all mixed signs, constraints and rotation identities | VERIFIED |
| Nonzero additive temporal frequency under E1 | OBSTRUCTION IDENTIFIED |
| Filled-axis winding under E2 | OBSTRUCTION IDENTIFIED |
| Source-free rotation vanishes under E3 boundary hypotheses | VERIFIED / OBSTRUCTION IDENTIFIED |
| Application of the global hypotheses to an arbitrary proposed geometry | CONDITIONAL |
| Complete circular stationary soliton under E4 | OBSTRUCTION IDENTIFIED |
| Noncircular or non-inheriting stationary-stress sectors | NOT YET ESTABLISHED |
| General stationary theorem, nonlinear stability, novelty, PRD readiness | NOT YET ESTABLISHED |

Observed execution results are recorded below only after running the commands.

Observed verification, 2026-11-10 UTC (Python 3, SymPy 1.14.0):

* Stage E: exit 0, **61 exact local checks PASS**, no FAIL. Includes direct
  four-dimensional curvature, both currents, full off-shell stress divergence,
  every spatial ADM component, rotation boundary product rule and angular
  momentum tail normalization. This is one investigator using independent
  computational routes, not an external verification team.
* Unchanged Stage D: exit 0, 27 exact checks PASS; unchanged Stage C: exit 0,
  all 23 exact identity checks PASS. Stage B, Stage A½ symbolic and Stage A½
  bounded: each exit 0 and final PASS. These retain their earlier scopes.
* Frozen-source provenance, `--package-only`: exit 0, 49 package hashes,
  20 archive members and pinned source copies PASS. Its historical recovery
  gate remains INCOMPLETE and physical tests remain NOT RUN, as before.
  No such tests are authorized or attempted in Stage E.
* `git diff --check`: exit 0. Baseline comparison against the exact Stage D
  commit, excluding only the three named Stage E paths: empty, exit 0.
  `git status --short` shows exactly those three new files. Thus all previous
  tracked files, including the original action, frozen evidence and every
  tracked Run 43 artifact, are unchanged; provenance package hashing supplies
  a separate content check for that package.

Only the three Stage E deliverables are to be committed and published on the
requested branch. Remote commit/tree verification is reported in the delivery
message (a commit cannot contain its own SHA). Stop after Stage E; further
research awaits separate authorization.
