# Stage H: anisotropic stationary spectral analysis

2026-11-10 UTC. One Codex agent. Branch
`research/stage-h-general-stationary-spectral`; exact Stage G parent
`ca30828547b5989c751b5c9639796377cdc1e3c9`.
Baseline tree `3c7f9c73bb4a3134ba84840452cde5020c4c0e28`, 1,123 tracked paths.
The companion `DEE_STAGE_H_GENERAL_STATIONARY_ASYMPTOTICS.md` defines
H1–H18, the unchanged theory, conditional end constraints and mass identity.

**Outcome C — Necessary conditions and residual spectral class.** The exact
anisotropic quadratic pencil and its long-range terms are derived. A general
absence-of-embedded-modes theorem applicable to the minimally specified
coupled end is **NOT YET ESTABLISHED**. No nonzero solution, counterexample
to nonexistence, or new anisotropic exclusion theorem is claimed. Stage G's
spherical-exterior obstruction remains valid within its original assumptions.
Failure to complete its generalization is a limitation of this proof, not an
existence result and not Outcome D. Neither Outcome A nor Outcome B is claimed.

## Exact spectral problem and functional class

Fix Omega real and nonzero, k=|Omega|. On a sufficiently far exterior end
E={R>R0}, h->delta,N->1,beta->0,f->f_v>0 make A uniformly positive definite
and bounded, w bounded above and below. Minimal asymptotic limits alone
need not control coefficient derivatives. For the following operator-domain
construction also assume the controlled smooth coefficient bounds in H6;
these bounds are stated hypotheses. The exact density equation is

    partial_i(A^ij partial_j z)
       +i Omega(2B^i partial_i+partial_i B^i)z
       +Omega² C0 z=0,
    A=w f(h^-1-beta beta/N²)/f_v,
    B=w f beta/(f_v N²), C0=w f/(f_v N²).                (S1)

Real symmetric A, real B,C0 and a complex profile are essential. The finite
matter energy argument in the companion report **derives** z in H¹(E) and
L²(E) with Euclidean or equivalent metric volume. The weak equation acts
H¹_loc to H^-1_loc; no radiation condition is implicitly part of it. Set

    H_Omega=-partial_i(A^ij partial_j)
                -i Omega(2B^i partial_i+partial_i B^i)
                +Omega²(1-C0),
    H_Omega z=Omega² z.                                 (S2)

For fixed Omega this is a differential operator in ordinary flat L²(dx).
Its coefficients depend on Omega as displayed, and through a coupled
solution also on the fields and metric. It is not one common linear
Hamiltonian eigenvalue problem obtained by varying Omega with the background
held on shell. This observation does not prevent applying a suitable theorem
to S2 separately at every fixed nonzero Omega.

A precise available realization on the artificial exterior with **Dirichlet**
boundary is the closed, lower-semibounded quadratic form on H¹_0(E)

    q_Omega[u]=integral_E [A^ij partial_i u* partial_j u
            +2 Omega Im(u* B.grad u)+Omega²(1-C0)|u|²]dx. (S3)

Uniform ellipticity, bounded B and bounded real 1-C0 give a lower bound by
Cauchy–Schwarz and the epsilon inequality; adding a sufficiently large
constant makes the form norm equivalent to H¹_0. Polarization defines a
Hermitian form. The representation theorem gives its Friedrichs self-adjoint
operator, with weak domain

    {u in H¹_0(E): H_Omega u in L²(E) distributionally}.

This is an actual self-adjoint realization under those assumptions; it
is distinct from formal symmetry of L on L²(w dx). The original z generally
has **nonzero trace at R0**. It must not be put in this Dirichlet domain
without proof. An exterior Rellich conclusion needed here must cover weak
solutions with arbitrary finite trace on the artificial interface, or be
applied after a justified global realization. Cutting off z introduces a
compactly supported source. Absence of eigenvalues for the Dirichlet exterior
by itself would not prove exterior vanishing of that source problem.
There is no physical reflecting surface at R0.

On a globally strictly timelike smooth complete slice with appropriate
uniform end bounds, compact interior coefficient bounds also permit a
semibounded form realization without the artificial boundary. If k becomes
spacelike, the principal form is indefinite and this Friedrichs elliptic
argument fails globally. The asymptotic end itself remains elliptic.

## Which established methods have and have not been applied

A bounded examination of relevant methods was made. No imported anisotropic
spectral theorem is used as a proven basis here without a verified matching
statement. This is a deliberate boundary on the conclusion, not a claim that
no such theorem exists in the mathematical literature.

* Flat Rellich arguments apply to an exact Helmholtz equation outside a
  compact set with the required decay/radiation integrals. S1 has tensor
  second-order and long-range frequency terms and is not that equation.
* Scalar Schrodinger criteria on a multiplicative potential do not alone
  control `partial_i((A-I)^ij partial_j)`. The shift is a first-order
  Hermitian differential term with a matched imaginary divergence, not
  an arbitrary complex absorbing potential.
* Long-range scattering and radiation estimates are a plausible route, but
  their differentiability, nonoscillation, weighted-space, possible trapping,
  solution-domain and boundary assumptions must be checked for S1. One must
  obtain an exterior Rellich uniqueness statement, not only a construction
  of outgoing solutions for a source.
* A positive commutator at infinity can imply estimates away from thresholds
  modulo compact errors. A Mourre estimate with a compact remainder alone
  does not exclude a finite set of embedded eigenvalues. Any eigenvalue
  exclusion needs an additional estimate/continuation argument. Interior
  trapping can obstruct uniform resolvent bounds without automatically
  producing a real L² eigenmode.
* Weighted elliptic Carleman estimates may upgrade decay and then give
  absence of embedded spectrum. They require control of principal jets and
  lower-order terms; no such estimate for the general minimally constrained
  Einstein–DEE end is proved by the symbolic verifier.

A verified bibliographic lead is H. Koch and D. Tataru, *Carleman estimates
and absence of embedded eigenvalues*,
https://arxiv.org/abs/math-ph/0508052. The arXiv abstract was retrieved in
this Stage H session. It states absence of embedded eigenvalues for a
Schrodinger potential in L^((n+1)/2) and extensions to variable coefficients,
long-range potentials and gradient potentials. The abstract alone does not
state the hypotheses of those extensions. A bounded attempt to retrieve
`https://arxiv.org/pdf/math-ph/0508052` returned HTML rather than a PDF;
`pdftotext` exited 1 and the combined extraction command exited 2. No full
extension theorem, equation number, or applicability certificate is invented.
The extension is an **unverified theorem lead**, not a proof ingredient.
In n=3 an unseparated nonzero Coulomb potential is not in L²(E), since
integral R² R^-2 dR diverges. Thus even the verified abstract's basic
integrability statement cannot directly handle that term. No conclusion
about the full paper's actual extended theorem follows from this observation.
The source download is unrelated to and not counted as a failed regression.

The needed exterior statement can be formulated independently of citations:
for the real coefficient class in S1/H6 (or a smaller justified class),
prove that every z in H¹(E) intersect L²(E) solving S1 has a zero exterior
germ, with arbitrary interface trace. This is the **missing spectral
statement**, not an assumption masquerading as a theorem. H6 is a concrete
candidate class with nonoscillatory differentiated 1/R leading jets;
it is not claimed here that its estimates are either optimal or known
insufficient. Still less is it claimed that finite ADM mass supplies H6.

## Long range, Coulomb phases, angular derivatives and amplitude

In H6 the full expanded equation is H8: A-I has an anisotropic O(1/R)
second-order piece, C0-1 an O(1/R) scalar piece, B=O(R^-p0), and amplitude
F=O(R^-s). The remaining A,C0 coefficients are O_k(R^-q0), q0>1, with
controlled differentiated decay. The symbol of the long-range problem
contains both `A^ij xi_i xi_j` and `2 Omega B^i xi_i`; it is not described
by the scalar tail alone. Angular derivatives of n1,h1 and b1 enter divergence
terms. Treating angular coefficients as radial constants loses those terms.

In the separately normalized isotropic mass gauge H13, set c=4m k².
For an exactly spherical exterior the harmonic ODE is

    y_lm''+[k²+c/R+r_l(R)]y_lm=0,
    r_l=O(R^-min(2,1+delta)) in L¹(dR),
    y_lm=R sqrt(N Hrad f) z_lm.                         (S4)

Hrad is the conformal radial factor; Hrad=1+m/R+... . Stage G assumed
O_2(R^-2) remainders, which is a special case. For fixed l the centrifugal
term is integrable. The comparison waves

    exp(+-i theta), theta=kR+c log R/(2k)

have exact residuals `[-c²/(4k²R²) -+ i c/(2kR²)]` times the wave, checked
again in the verifier. Integrable ODE variation-of-constants then gives
nonzero limiting oscillation amplitudes for every nonzero solution. Their
mean squared magnitude is positive and contradicts L²(dR). This is the
existing controlled spherical obstruction; all harmonics are treated one
at a time, not in a uniform large-l approximation. An arbitrary tensor
perturbation couples harmonics and replaces that scalar ODE with an
unbounded angular operator. The scalar comparison fundamental-matrix proof
cannot simply be reused with a finite-dimensional Gronwall bound.

A direct phase transformation also shows the limitation. For a flat
Coulomb comparison equation, put z=exp(i phi)u. Then exactly

    exp(-i phi)(Delta+k²+c/R)exp(i phi)u
      =Delta u+2i grad phi.grad u
            +[k²+c/R+i Delta phi-|grad phi|²]u.          (S5)

The unit-modulus choice phi=c log R/(2k) has |grad phi|²=O(R^-2),
so it **does not remove c/R from the full PDE**. It corrects the eikonal
phase of a selected outgoing or incoming wave when combined with its
leading exp(+-ikR), and those two channels require opposite signs. In a
variable metric the transformed principal tensor A remains A and there are
additional A-grad-phi and B-grad-phi terms. A single scalar phase cannot
eliminate a general anisotropic principal perturbation or all propagation
directions. Modified scattering instead solves the appropriate directional
eikonal and transport equations and needs estimates for its remainders.
No such global construction is assumed by saying 'remove the Coulomb term'.

With the H6 amplitude bounds, f variation supplies the explicitly displayed
shorter-range F terms and their derivatives, rather than a new phase mass.
Without those bounds H5 has only the energy-controlled bilinear source;
the amplitude can then be another missing source of coefficient-jet control.
The spin term in the conditional rest gauge is
`+4G Omega J.L/R³`, not a scalar central attraction. For general beta it
has the first-order form in S1; its sign depends on profile and frequency.

For a positive-commutator diagnostic, let D=-i(x.grad+3/2) and
V_Omega=Omega²(1-C0). The formal dilation commutator has principal density
`2A-x.grad A`, potential `-x.grad V_Omega`, and shift transport coefficient
`B-x.grad B` in the same Hermitian first-order combination. The verifier
checks the full local one-dimensional product identity, including shift,
as a bounded algebraic diagnostic; no full PDE estimate is claimed from that
test. For A=I and attractive V_Omega=-c/R, the potential contribution is
V_Omega, which is negative for c>0. Even when the principal end term dominates
far away, compact interior errors prevent a claimed global coercivity identity
without further work. Positive matter energy is not this commutator positivity.

## Radiation and weights

Use L²_sigma(E)={u: <R>^sigma u in L²(dx)} and H¹_sigma with the same
weight on first derivatives. Finite proper matter energy supplies unweighted
H¹ and L², not arbitrary positive weights or a pointwise Sommerfeld condition.
For the normalized Coulomb channels a candidate modified outgoing condition
controls `(partial_R-i[k+c/(2kR)])(R z)` in appropriate angular/weighted
norms. The incoming channel has the opposite sign. Its precise norm and
required geometry must be stated by a scattering theorem. An L² hypothesis
can imply weak averaged radiation smallness, but does not itself prove the
asymptotic channel expansion or completeness of that expansion for S1.

Multiplying S1 by z* and taking the real part with exhaustion cutoffs gives
the exact no-source balance

    integral N f(|Dz|²-|p|²)dvol_h=0                   (S6)

when the integrals converge and the cutoff errors vanish. In the normalized
end those errors are bounded by tail L² norms and O(1/R) cutoff gradients.
The imaginary Green flux is conserved between interfaces; an L² mode has
vanishing flux along a suitable exhaustion. Zero net flux also allows a
standing combination of incoming and outgoing channels; it cannot replace
the positive mean-amplitude argument in S4. S6 equates two nonnegative
quantities but does not set either to zero. Shift contributions inside p
are crucial. A positive-flux argument needs an outgoing condition or a proven
radiation expansion, neither silently added here.

Embedded eigenvalues are the precise concern: Omega²>0 lies in the massless
asymptotic continuum. Attractive Coulomb negative-energy bound levels of a
fixed Hamiltonian do not address that question; moreover V_Omega depends
on Omega. Known oscillatory long-range counterexamples in other scalar
spectral classes motivate derivative hypotheses, but no such potential is
constructed here, and none is claimed to satisfy the Einstein–DEE equations.
A constructed test profile or an artificial external potential would not
resolve the coupled question.

## Coupling, additional assumptions and continuation

H9 supplies necessary constraints; H10–H11 prove useful restrictions on leading
angular data **after** H6 is assumed. H6 plus finite matter energy and DEC
force the R^-3 stress jet to vanish, hence n1 is constant and h1 obeys the
linearized static-vacuum tensor equation. This is stronger than treating all
angular functions as free. It does not derive H6, an admissible global
isotropic/rest gauge, its higher differentiated remainders, a spin tail,
or the missing exterior spectral estimate. The lapse source is U-V and
includes beta.DK; it has no fixed sign. H18 is a valid noncircular charge/mass
balance under its boundary conditions and gives no contradiction for nonzero
Omega. The positive-mass theorem constrains ADM mass, not absence of embedded
spectrum for S2.

The weakest assumptions sufficient for a new proof have **not** been
established, and no optimal decay threshold is asserted. A concrete program
that would close the gap must separately supply:

1. An admissible stationary end and weighted regularity theorem from the
   coupled equations, or explicit extra assumptions such as H6 and a
   justified rest normalization H13 with differentiated remainders.
2. A proved exterior Rellich/Carleman or radiation theorem for exactly that
   real tensor/shift long-range class and unweighted H¹ intersect L² weak
   solutions, covering arbitrary artificial-interface trace. A theorem's
   function-space hypotheses must follow from the energy or an additional
   justified decay upgrade.
3. A valid inward continuation argument for the actual interior geometry.

These are distinct missing steps, not a circular theorem saying that modes
are absent if modes are absent. H6 gives explicit coefficient estimates for
step 2; it is recorded as a target class for a theorem, not as a proved
sufficient theorem. An already rigorous sufficient geometric class is the
exactly spherical exterior of Stage G with its globally timelike interior.
Stationary admissible coordinate pullbacks of that class can have anisotropic
coefficients: transport the profile by H12 and transport the energy by the
bounded coordinate/slicing equivalence. This is an immediate coordinate
invariance corollary of Stage G, not a new physically general anisotropic
obstruction and not a reason to classify Stage H as Outcome B.

If an exterior theorem eventually proves z=0 on an open region, on any
connected region with fP real symmetric locally uniformly positive definite,
principal coefficients locally Lipschitz, and the complex first/zero-order
coefficients locally bounded, standard local weak elliptic unique continuation
is the applicable input. N and f must remain positive with positive local
lower bounds, h nondegenerate, and all these bounds must hold on the compact
balls used for continuation. Smooth strictly timelike stationary coefficients
satisfy this; uniform ellipticity across all of infinity is unnecessary for
this local step. On a real/imaginary two-component representation the principal
operator is the same scalar operator on both components; the lower-order
matrix couplings obey a bound C(|grad u|+|u|). Scalar Carleman inequalities
can be summed over the two components and absorb these bounded couplings.
The local continuation conclusion propagates exterior vanishing through a
chain of overlapping coordinate balls on the connected timelike component.
This is the standard conditional elliptic input used by Stage G, not a new
symbolic or independently peer-reviewed PDE theorem. Stage G cites
N. Aronszajn, J. Math. Pures Appl. 36 (1957), 235–249. Its original paper
was not retrieved in Stage H; the inherited citation is not represented as
a newly verified external theorem statement. The smooth elliptic continuation
step is recorded with its hypotheses and cannot repair the missing exterior
spectral proof.

At an ergosurface the elliptic principal form degenerates; inside, it can
have mixed type. Thus even conditional exterior vanishing can currently be
continued only in the connected strictly timelike region. It does not yield
a zero profile throughout an arbitrary ergoregion. The full original field
system is hyperbolic in spacetime, but vanishing on an exterior worldtube is
not full zero Cauchy data. Hyperbolic finite-propagation/Cauchy uniqueness
alone therefore does not fill this gap. Possible continuation across a
noncharacteristic surface needs a suitable pseudoconvexity/analyticity or
other rigorous theorem; none follows here from smoothness, stationarity or
DEC. No additional analyticity is assumed. The ergoregion sector remains
**NOT YET ESTABLISHED**. This is neither a solution nor a counterexample.

## Scientific classification and publication limits

| Claim | Status |
|---|---|
| Exact anisotropic stationary operator and formal Green identity | VERIFIED |
| Positive matter energy implies H¹ intersect L² on the normalized end | VERIFIED under the stated end bounds |
| Differentiated H6 expansion from finite ADM mass alone | NOT YET ESTABLISHED |
| H7–H8 jets and leading Einstein restrictions after H6 | VERIFIED / CONDITIONAL |
| Rest-gauge scalar Coulomb coefficient -4G M_ADM Omega²/R | CONDITIONAL; exact first-jet algebra VERIFIED |
| Same scalar coefficient in every stationary coordinate gauge | OBSTRUCTION IDENTIFIED: principal tensor also changes |
| Noncircular Komar cancellation H16 | VERIFIED |
| ADM identity without boundaries and mass-equality conditions | OBSTRUCTION IDENTIFIED; H18 is the corrected general form |
| New broad anisotropic positive-frequency exclusion | NOT YET ESTABLISHED |
| Interior continuation through a general ergoregion | NOT YET ESTABLISHED |
| General stationary Einstein–DEE mode or universal no-go | NOT YET ESTABLISHED |

**Strongest defensible mathematical conclusion:** the original phase equation
requires a positive-frequency L² mode of the explicitly real, formally
Hermitian variable-coefficient pencil S1. The leading stationary Einstein
constraints restrict any controlled 1/R geometry, but no derivation from the
minimal global assumptions to a fully applicable anisotropic exterior
absence theorem is completed here. The noncircular Komar identity is valid
with its boundary/end qualifications and supplies no missing coercivity.
Stage G remains an obstruction within its spherical-exterior, elliptic class.

**Single most important unresolved question:** can the coupled finite-energy
stationary Einstein–DEE equations be shown to place their exterior density
pencil in a rigorously verified absence-of-positive-L²-modes class, with the
necessary differentiated decay and admissible gauge? Even a positive answer
would leave interior ergoregion continuation to address for a universal result.

**Conservative Physical Review D assessment:** not ready as a claimed general
stationary exclusion or existence result. The exact reduction, conditional
constraint analysis and mass reconciliation are useful research records;
they do not establish novelty, a complete global spectral proof, stability,
independent mathematical review or physical realization. Internal exact
symbolic checks are not independent peer review. No PRD manuscript or
submission, simulations, theory modifications or Stage I is authorized here.

## Verification scope and actual execution records

The standalone `DEE_STAGE_H_VERIFIER.py` uses exact SymPy arithmetic, no
random sampling, no floating tests, no network, and no file outputs. The
named-check count counts identities, some with several tensor components.
It verifies general ADM inversion/determinant, the exact reduced pencil,
principal form, anisotropic Green flux, density first variations, metric/
lapse/frequency/shift/amplitude terms, spatial gauge strain and linearized
Ricci constraints, time-profile gauge phase, exact spherical Stage G reduction,
Coulomb comparison residuals, original-stress Komar/Noether cancellation,
boundary product rules and a documented local dilation commutator.
It does not certify global PDE theorems, optimal coefficient assumptions,
Einstein solution existence, ADM–Komar equality without hypotheses, or
continuation in mixed type. H10's use of integrated DEC and the sphere
harmonic argument is analytical, not part of the symbolic counter.

Exact algebraic checks, analytical deductions with stated hypotheses,
standard mathematical inputs, and unverified literature leads are
separated above. Full actual verification output and independent exact-tree
preservation results follow; temporary execution records stay in /tmp.

### Final observed verification runs, 2026-11-10 UTC

Python 3.12.14, SymPy 1.14.0. These are fresh executions, not historical counts.
Stage H: **43** named exact local checks; G: **60**; F: **76**; E: **61**; D: **27**; C: **23** named PASS identities. Stage B and both Stage A½ verifiers retain their unnumbered final PASS. All ten final verification commands exited **0**. The record below is the complete stdout and stderr for the final required runs; empty stderr is stated explicitly.

Command: `python DEE_STAGE_H_VERIFIER.py`; exit code **0**.

Stdout:

```text
SymPy 1.14.0
PASS: general ADM inverse
PASS: ADM determinant by unit shear
PASS: ADM diagonal-block determinant
PASS: principal quadratic form
PASS: orthonormal principal eigenvalues
PASS: general stationary covariant phase reduction
PASS: density pencil identity
PASS: full anisotropic Green flux
PASS: shift skew transport product rule
PASS: first variation inverse spatial metric
PASS: first variation spatial determinant
PASS: first variation volume sqrt determinant
PASS: first variation spacetime density
PASS: first variation spatial density
PASS: first variation frequency density
PASS: first variation shift density
PASS: first variation shift squared principal term
PASS: metric differential versus scalar correction
PASS: amplitude stiffness exact expansion
PASS: exact division removes multiplicative stiffness
PASS: isotropic mass spatial density cancellation
PASS: isotropic mass Coulomb frequency
PASS: areal gauge trace
PASS: areal gauge radial principal correction
PASS: areal versus isotropic spatial gauge
PASS: linearized Ricci pure spatial gauge
PASS: linearized stationary mass equation
PASS: linearized momentum operator
PASS: stationary time gauge shift is momentum pure gauge
PASS: time gauge covariant profile derivative
PASS: exact stationary time coordinate inverse metric
PASS: full PDE phase conjugation retains Coulomb term
PASS: exact spherical Stage G operator
PASS: exact radial density conjugation
PASS: Coulomb logarithmic phase residual 1
PASS: Coulomb logarithmic phase residual -1
PASS: general shift Komar charge cancellation
PASS: trace reversed Komar integrand
PASS: normal Noether density from contravariant time current
PASS: momentum contracted shift boundary product
PASS: Komar surface from two ADM integrations
PASS: phase real energy boundary product
PASS: local dilation commutator including shift
PASS: Stage H 43 exact local checks; analytical and unresolved claims are not symbolic tests.
```

Stderr: empty.

Command: `python DEE_STAGE_G_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python DEE_STAGE_F_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python DEE_STAGE_E_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python DEE_STAGE_D_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python DEE_STAGE_C_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python DEE_STAGE_B_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python DEE_STAGE_A_HALF_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python DEE_STAGE_A_HALF_BOUNDED_VERIFIER.py`; exit code **0**.

Stdout:

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

Stderr: empty.

Command: `python evidence/einstein_dee_sources/verify_evidence.py --package-only`; exit code **0**.

Stdout:

```text
Git preservation NOT CHECKED in package-only mode
PASS: 49 package hashes, 20 archive members, pinned source copies
Recovery gate: INCOMPLETE. Physical tests: NOT RUN; approval required.
```

Stderr: empty.

Frozen package: **49 hashes, 20 archive members and pinned copies PASS**, exit 0. Recovery gate remains **INCOMPLETE**; physical tests **NOT RUN**. The printed approval note is historical verifier output, not a request to run new physical tests. Package-only mode does not verify Git preservation.

### Independent exact-parent preservation

Command: `python /tmp/dee-stage-h/preserve.py`; exit **0**. Complete stdout:

```json
{
  "baseline": "ca30828547b5989c751b5c9639796377cdc1e3c9",
  "baseline_tree": "3c7f9c73bb4a3134ba84840452cde5020c4c0e28",
  "head": "ca30828547b5989c751b5c9639796377cdc1e3c9",
  "tracked_paths_checked": 1123,
  "run43_named_paths": 233,
  "content_or_mode_mismatches": [],
  "index_mismatches": [],
  "head_baseline_mismatches": [],
  "additions": [
    "DEE_STAGE_H_ANISOTROPIC_SPECTRAL_ANALYSIS.md",
    "DEE_STAGE_H_GENERAL_STATIONARY_ASYMPTOTICS.md",
    "DEE_STAGE_H_VERIFIER.py"
  ],
  "exactly_three_authorized_additions": true,
  "unauthorized_index_paths": [],
  "unauthorized_head_paths": []
}
```

Stderr: empty.

The comparison enumerates every recursive parent entry and recomputes each worktree Git blob ID from exact file bytes, or symlink-target bytes. It checks Git executable/symlink modes, the index, and HEAD against the exact baseline. It also walks all worktree files (including ignored files, excluding only .git) and requires exactly the three authorized additions. The recorded complete NUL-delimited baseline tree is compared again to `git ls-tree -rz` at the immutable parent. All **1,123** historical tracked paths, including **233 Run 43 named paths**, match; no historical file, mode or symlink target changed. The full-tree check is independent of package hashes.

`git diff --check`: exit **0**, stdout and stderr empty before staging. Final staging checks and the same exhaustive comparison are repeated before commit. Only these three files are authorized:

- `DEE_STAGE_H_GENERAL_STATIONARY_ASYMPTOTICS.md`
- `DEE_STAGE_H_ANISOTROPIC_SPECTRAL_ANALYSIS.md`
- `DEE_STAGE_H_VERIFIER.py`

The requested branch is committed only after these checks succeed. Remote branch SHA, commit parent, all three remote blob IDs, unchanged historical tree entries, clean index/worktree and unrelated local refs are checked after push and reported in the delivery message. A commit cannot contain its own SHA; no postcommit report edit is required. No main merge, force push, amend, cleanup, history rewrite, Stage I, simulation or manuscript work is performed. **Stop after Stage H publication.**
