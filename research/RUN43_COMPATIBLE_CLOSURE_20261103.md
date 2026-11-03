# Run 43: executed compatible-closure alternatives, November 3

**Objective not met. No candidate is promoted.** This is new implementation and
execution, not a repetition of the prior saved-field studies. All four attempts
completed n=8,16,32,64 with original matter. The physically justified direct
closures fail the unchanged frozen gate. The initial angular metric scaling
error is preserved explicitly, rather than counted as a valid improvement.

Branch: `development/run43-compatible-hyperviscosity-20261103`, based on verified
`6819962` and its parent `f110d62`. No merge or journal submission was made.
The November 2 result was not used as replacement physical validation.

## Continuum equations and complete boundary construction

The unchanged constraints are

    Delta psi + |LW|^2 psi^-7/8 + 2 pi G E(psi) psi^5 = 0
    div(LW) = 8 pi G psi^6 S
    E(psi) = temporal/2 + spatial psi^-4/2 + U
    LW_ij = partial_i W_j + partial_j W_i - 2 delta_ij div W/3.

All Run38 fields, normal velocities, kinetic mixing, G=.001 and the inverse
physical metric in E are retained. Dirichlet psi=1,W=0 at .05,.25 and periodic
angles are exact. No residual is absorbed into E or S. The independent 4D
stress projections are checked on each new solution.

For Cartesian scalar component f, use analytic toroidal-frame Hessians
(q=1+r cos theta):

    Hrr = f_rr
    Htt = f_tt/r^2 + f_r/r
    Hpp = f_pp/q^2 + cos(theta) f_r/q - sin(theta) f_t/(r q)
    Hrt = f_rt/r - f_t/r^2
    Hrp = f_rp/q - cos(theta) f_p/q^2
    Htp = f_tp/(r q) + sin(theta) f_p/q^2.

Then Hess_cart(f)=sum_a ea ea^T Haa + sum_(a<b)(ea eb^T+eb ea^T)Hab.
For each W component, Delta_L W_i=trace(Hess W_i)+sum_j Hess(W_j)_ij/3.
Thus radial, angular, mixed and connection terms are assembled together;
there is no fitted decomposition or cancellation subtraction. In particular
radial second derivatives are DIRECT, never compositions of endpoint D1.
The radial boundary closures are polynomial moment systems, with weights w
satisfying sum_j w_j (j-i)^k = k! delta_(k,d)/h^d. Six radial D2 points
in the fourth-order assembly reproduce degrees through five. The standard
second-order variant uses four D2 points, reproducing degrees through three.
Both remove the O(h) first-interior-row composition term from the SOLVER.
They do not change the frozen audit's derivative matrices.

The coordinate-KO candidate adds

    -sum_a h_a^2/(4 a_a^2) (I+ea ea^T/3) D_aaaa W,
    (a_r,a_theta,a_phi)=(1,r,q).

The coefficient 1/4 follows the centered interior identity
D2-D0 D0=-h^2 D4/4. Its nonpositive constant-coefficient symbol damps the
Nyquist mode. Six-point radial D4 closures reproduce degrees through five,
so the added term is uniformly O(h^2) on smooth fields, including boundary
rows, without extra physical boundary data. Angular physical spacings are
r*h_theta and q*h_phi. The initial implementation mistakenly used a_a^-4
instead of a_a^-2; `run43_compatible_closure_metric4.py` and its entire execution
are retained. It is explicitly rejected as the intended mesh-scaled method.

The no-dissipation second-Hessian candidate is the standard second-order direct
operator with all six Hessian terms. The fourth alternative uses the original
tensor-divergence operator plus h_r^4/16 (I+er er^T/3) D6 W. Its fixed coefficient
matches angular-filter Nyquist damping; eight-point radial D6 closures need no
additional boundary condition. This last alternative retains the original O(h)
composition defect and tests stability only. It is not a uniform boundary
accuracy remedy. No parameter was fitted to an order or a field.

The Hamiltonian retains the established conservative scalar face flux and its
fixed .25 angular fourth difference; its LW uses each candidate's gradients.
This is a consistent numerical operator, not a change to the physical sources.
The same nonlinear coupling, stopping conditions and Krylov tolerances apply.
All actual vector residuals independently satisfy their stated tolerances.

## All twelve orders from a separate implementation

The original frozen audit blob remains
`5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`. Every one of the twelve orders must
be finite and >=1.5. Band inclusion, interior-row selection, volume weights
r*q*psi^6 and H/M2 definitions are unchanged. `audit_run43_arrays.py` reloads
and hashes arrays in a separate process. The new
`verify_run43_closure_fields.py` independently reconstructs H, M2 and stress
projections using explicit array differences; it does not import Audit or
solver derivative/tensor code. Its norms agree with the frozen audit and
produce the orders below. Agreement is an implementation cross-check, not
external scientific replication.

### compatible-closure-20261103: 10/12, rejected

| Band | Constraint | 8->16 | 16->32 | 32->64 |
|---|---|---:|---:|---:|
| inner | H | 1.740848297 | 1.929827016 | 1.996176044 |
| inner | M2 | 0.252079856 | 1.132307944 | 1.672118109 |
| outer | H | 1.745664711 | 1.963989447 | 2.018377652 |
| outer | M2 | 1.528303103 | 1.755251217 | 1.919946436 |

### coordinate-ko-20261103: 9/12, rejected

| Band | Constraint | 8->16 | 16->32 | 32->64 |
|---|---|---:|---:|---:|
| inner | H | 1.740849026 | 1.929825857 | 1.996175145 |
| inner | M2 | 0.392887689 | 1.148509368 | 1.354471263 |
| outer | H | 1.745667846 | 1.963990187 | 2.018376424 |
| outer | M2 | 1.902861480 | 1.981727586 | 1.996286747 |

### second-hessian-20261103: 9/12, rejected

| Band | Constraint | 8->16 | 16->32 | 32->64 |
|---|---|---:|---:|---:|
| inner | H | 1.740847734 | 1.929825509 | 1.996175316 |
| inner | M2 | 0.508241329 | 0.852441285 | 1.092221118 |
| outer | H | 1.745664036 | 1.963989187 | 2.018377094 |
| outer | M2 | 1.844384409 | 1.967535936 | 1.991607324 |

### radial-filter-20261103: 9/12, rejected

| Band | Constraint | 8->16 | 16->32 | 32->64 |
|---|---|---:|---:|---:|
| inner | H | 1.740847607 | 1.929825592 | 1.996175340 |
| inner | M2 | -0.013997314 | 0.793939001 | 1.224283403 |
| outer | H | 1.745663884 | 1.963989080 | 2.018376564 |
| outer | M2 | 1.539939562 | 1.500286568 | 1.514581409 |

## Correction-specific failure evidence and limits

The direct closures pass independent analytic radial tests, including the first
interior rows; second-order candidate tests also include angular and mixed
contributions using a Cartesian polynomial with known vector Laplacian.
The coordinate-KO radial first-row finest tested orders are 2.087 inner and
1.865 outer. Small-grid complete eigenvalue checks find no zero or growing
mode for either direct candidate. All 43 regression tests pass. These finite
checks are not a proof of uniform stability on arbitrary grids.

On the NEW second-Hessian n=64 field, supplementary original-source inner M
is 8.5644144e-6, while frozen inner M2 is 1.9897416e-4. The exact frozen-minus-
matched endpoint term has norm 2.0061843e-4. Thus the large audit/closure
mismatch persists on the newly corrected field. Norms are not additive: no
cancellation is fitted or removed. On the new coordinate-KO field these three
norms are 7.6813408e-5, 2.2695556e-4 and 2.0394397e-4. Dissipation adds material
coarse-grid error without curing the frozen endpoint mismatch. The radial-sixth
filter is also pre-asymptotic: its first inner-M2 order is negative, while its
new n=64 supplementary inner M is 7.6290051e-4. It does not repair the closure.
`endpoint-discrepancy.json` in each candidate records the exact endpoint term
on all four new fields; this is a specific validation of these corrections,
not a repeat of the old full factorial budget.

None of these facts is an impossibility proof. Nor does a successful test of
operator accuracy warrant relabeling any failed gate as a pass. The prior
fourth-vector candidate remains more accurate under the supplementary metric
than these alternatives; no replacement is claimed. Fixed-source conformal
initial constraints on a finite shell do not constitute full Einstein-DEE
matter evolution or physical DEE validation.

## Preservation and reproducibility

`results/compatible-closure-20261103-evidence/` contains all execution and audit
logs (including initial checker errors and their corrected runs), summary of
all checks, environment and provenance, true-residual verification and complete
SHA256SUMS. Each of the four candidate directories contains raw fields, every
completed nonlinear checkpoint, full GMRES histories, executed source snapshots,
original-source audits, frozen-code audits and independent implementation audits.
All 360 pre-existing files were hashed before and after; zero changed. All
105 prior NPZ archives are unchanged. Original solvers and frozen audit were
not edited. The development branch contains new files only.

Replay from repository root with NumPy/SciPy as recorded in environment.json:

    DEE_REPLAY_TAG=unique-new-tag bash research/reproduce_run43_compatible_closure.sh

The script refuses existing output directories and retains failing audit exit
codes while completing subsequent audits. Both frozen-code and independent
implementation gates return exit 1 for scientific failure, and write all twelve
checks. The source closure script names and command environment are explicit;
no November 2 field or modified source is loaded for a new solve.

No candidate is merged, promoted, or submitted for publication. An original-
source 12/12 correction remains unestablished by this execution.
