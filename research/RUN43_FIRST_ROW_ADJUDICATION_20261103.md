# Run 43: first-row consistency adjudication

Base: `7ae549a90c94a979132affff255cabd9471f7146`.
Branch: `research/run43-first-row-consistency-20261103`.

**Scientific blocker: no additional production correction is supported by this
evidence.** The frozen finite-grid inner-M2 test rejects exact continuum
manufactured solutions and much more accurate fourth-order numerical solutions,
while accepting momentum convergence of less accurate tensor solutions. This is
not a theorem that no consistent original-source scheme can pass all twelve
checks. No criterion changes and no candidate is promoted. The existing
original-source fourth-vector candidate remains **9/12 FAIL**.

## New work and scope

The previous endpoint identities, factorial budgets, source projections and
production all-twelve audits were not rerun. New work comprises a complete
first-row matrix certificate, manufactured **boundary-value solves with absolute
solution errors**, a radial energy proof, independent high-precision Cartesian
verification of the analytic forcing, and a row-local physical residual identity.
Earlier manufactured checks measured operator truncation without solving these
analytic boundary-value problems. No original solver or physical source changes.

## Exact first-row operator difference

Write G_i=sum_alpha a_i,alpha D_alpha, with coefficient columns
(er_i,et_i/r,ep_i/q), q=1+r cos(theta). The original expanded operator is

    E_i = L W_i + G_i sum_j G_j W_j/3.

The frozen audit applies the unexpanded tensor operator

    T_ij = delta_ij sum_k G_k G_k + G_j G_i - 2 G_i G_j/3.

The already-established identity is
(T-E)_i=(sum_j G_j G_j-L)W_i+sum_j [G_j,G_i]W_j. The experimental tensor solver
also adds its fixed angular dissipation; it is not identical to either finite-grid
operator. All represent the same continuum Delta W+grad(div W)/3.

The NEW certificate compares T with direct second-order Cartesian Hessian H,
using identical first derivatives. Define direct pure chart derivatives Q and

    B_alpha = D_alpha D_alpha - Q_alpha,alpha
    C_alpha(b) = D_alpha diag(b)-diag(b) D_alpha-diag(partial_alpha b).

Coefficient derivatives are analytic. The exact difference of G_i G_j from the
direct Cartesian Hessian is

    K_ij = sum_alpha diag(a_i,alpha a_j,alpha) B_alpha
         + sum_alpha,beta diag(a_i,alpha) C_alpha(a_j,beta) D_beta.

Mixed chart derivatives commute exactly as tensor-product matrices. Therefore

    (T-H)_ij = delta_ij sum_k K_kk + K_ji - 2 K_ij/3.

This retains every angular, radial, mixed, geometry and derivative-composition
term, including finite-grid products of derivative errors. No Taylor term or
fitted cancellation is substituted for the full matrices. The certificate
includes all input coefficients, including both boundary planes, before
Dirichlet elimination. At n=8 and 16 the matrices have shapes 192x1728 and
768x13056. Relative coefficient identity errors are **5.16e-16 and 5.94e-16**.
Complete radial, angular, Leibniz and difference sparse matrices are archived.

At row 1, independently of W,

    (Dr Dr f)_1 = (3 f0-5 f1+f2+f3)/(4h^2)
    (Dr Dr-Drr)f_1 = (-f0+3f1-3f2+f3)/(4h^2) = h delta3 f/4.

The radial vector term is exactly (I+er er^T/3)h delta3 W/4. Its coefficient
identity errors are 4.55e-13 and 2.73e-12, on coefficients of magnitude 1600
and 6400. Angular and Leibniz matrices are material and are not discarded.

## Continuum consistency and boundary behavior

For smooth Cartesian components, the first interior radial composition satisfies

    Dr Dr f_1 = f''(r1)+h f'''(r1)/4+5h^2 f''''(r1)/24+O(h^3).

At the inner endpoint, Dr Dr f_0=f''(r0)+3h f'''(r0)/4+O(h^2); outer endpoints
reverse odd-term signs. These operators are **consistent**, but have first-order
second-derivative truncation near the radial boundary. Angular composition and
mixed/product-rule terms are normally second order on smooth fields. Direct H2
removes the O(h) radial term; H4 has fourth-order boundary closure. Their
continuum vector constraints and physical sources are identical.

Dirichlet W=0 makes tangential derivatives vanish on the boundary, but does not
make the third normal derivative vanish. Neither the fixed original matter nor
the Dirichlet data supply that extra condition. Excluding the endpoint from the
norm does not remove its derivative from the first interior stencil.

One O(h) row contributes O(h^(3/2)) to a fixed-band volume RMS asymptotically.
This does not guarantee p>=1.5 on finite pairs. In the exact examples below,
f'''=6ab/r^4. The first row approaches a region of substantially larger third
derivative under refinement, changing the leading coefficient on the four
specified grids. Orders below 1.5 need neither instability nor incorrect forcing.

## Independent manufactured boundary-value solves

Set a=.05,b=.25,c=(.3,-.4,.5), f(r)=a+b-r-ab/r. Use exact fields
W=c f(r) and W=c f(r) sin(theta). Both have exact homogeneous radial Dirichlet
data and are smooth on the shell; the singularity lies outside it. The second
field has genuine mixed radial-angular derivatives.

Forcing is the analytic continuum Cartesian vector Laplacian, not a discrete
stencil applied to the exact field. Independently use
Hess r=et et^T/r+cos(theta) ep ep^T/q. For mixed F=z g(r), g=f/r,

    Hess F = g'(er ez^T+ez er^T)+z(g'' er er^T+g' Hess r)
    Delta_L(cF) = c trace(Hess F)+(Hess F)c/3.

A separate test differentiates F(x,y,z) directly using mpmath at 50 digits,
at six field/point combinations, agreeing within rtol=2e-12, atol=2e-11.
No production S, E or matter coefficients are modified. Manufactured forcing
is solely an operator test and is never substituted for physical acceptance.

**Sixteen new linear solves completed:** two fields, two existing operators,
n=8,16,32,64. All true residuals meet the test tolerances rtol=1e-11,
atol=1e-13. Raw archives contain exact fields, solved fields and analytic forcing;
JSON includes full GMRES histories, residuals, tolerances and hashes. The split

    T Wnum-forcing = (T Wexact-forcing)+T(Wnum-Wexact)

closes numerically and separates audit consistency error from solution effects.
Norms of these terms are not treated as additive scalars.

Inner-M2 orders on **exact continuum fields**:

| Field | 8→16 | 16→32 | 32→64 |
|---|---:|---:|---:|
| radial | 0.309947903 | 0.795480313 | 1.104801455 |
| mixed | 0.181239537 | 0.741036929 | 1.080163082 |

All six fail 1.5. This is a counterexample to a universal finite-grid accuracy
interpretation, not a replacement physical acceptance run.

At n=64, inner vector RMS values are:

| Field | Solver | Absolute solution error | Audit on solved | Audit on exact |
|---|---|---:|---:|---:|
| radial | tensor | 2.145608782e-4 | 1.730472387e-5 | 1.139112621 |
| radial | fourth | 6.786563754e-8 | 1.142792824 | 1.139112621 |
| mixed | tensor | 1.035585331e-4 | 1.080665660e-2 | 0.808032055 |
| mixed | fourth | 4.852232442e-8 | 0.810580261 | 0.808032055 |

The accurate fourth solution's audit residual approaches the exact field's audit
truncation error. Larger tensor solution errors have discrete derivatives that
counterbalance audit consistency error. This cancellation is measured, not
fitted. The mixed tensor solution's first absolute-error order is 1.436, while
its audit-residual orders are 1.907,1.984,1.988. Finite audit ratios do not reliably
rank absolute solution accuracy. Very high final fourth-solution orders in
this example reflect pre-asymptotic cancellation/superconvergence; no universal
order above four is claimed. All outer-band values and orders are in report.json.

## Stability: proof and limits

The O(h) composition defect does not itself imply instability. For
L=(Dr Dr)_II after Dirichlet elimination, let A=-4h^2 sym(L), and S be the
odd/even Dirichlet tridiagonal graph Laplacians, with edges at index distance 2.
For every even n>=8, A-S/2 is a diagonally dominant symmetric M-matrix with
nonnegative row sums, hence positive semidefinite. The longest chain has n/2
nodes, giving

    lambda_min(-sym L) >= sin(pi/(n+2))^2/(2h^2) >= 32, h=.2/n.

The last inequality follows from sin x>=2x/pi and n/(n+2)>=.8. This proves
uniform radial energy dissipation and an inverse Euclidean norm bound <=1/32.
The radial vector symbol I+er er^T/3 has eigenvalues 1,1,4/3 and er independent
of r, so the pure radial vector part inherits the bound. The executable graph
certificate and spectra are checked through n=256. Observed original radial
energy lower bounds rise from 212.50 to 245.83. Fourth-order radial energy
bounds are positive on those tested grids too.

This is **not** a uniform stability theorem for the full variable-coefficient
toroidal vector operator. Existing complete small-grid spectra remain valid;
new full-grid manufactured solves show convergent absolute errors with small
true residuals. This is tested evidence, not a proof for arbitrary grids or
sources. Neither radial checkerboard singularity nor Krylov breakdown explains
the measured failures.

## Original-source first-row residual identity

Only the unresolved row decomposition was evaluated on existing original-source
fourth-vector fields. No production solve or all-twelve audit was repeated.
With b=8 pi G psi^6 S and H the actual solved fourth-order operator,

    M_frozen = psi^-10(H W-b)+psi^-10(T W-H W).

The first term is physical algebraic solve residual; the second is operator
mismatch. Both use the same original physical sources and psi factors. The
identity closes to <=5.15e-18 across all interior rows. At the first row:

| n | Physical solve RMS | Operator-mismatch RMS | Frozen M2 RMS |
|---:|---:|---:|---:|
| 8 | 3.354811462e-13 | 1.252558428e-3 | 1.252558427e-3 |
| 16 | 5.685430976e-13 | 4.213922070e-4 | 4.213922072e-4 |
| 32 | 2.008737586e-14 | 6.119317347e-4 | 6.119317347e-4 |
| 64 | 1.032513354e-13 | 6.833769963e-4 | 6.833769963e-4 |

Thus the first-row frozen residual is explained by operator mismatch to
roundoff, not unresolved linear solve error. This identity alone does not prove
zero continuum solution error. The manufactured solves establish how an accurate
solution retains the mismatch; preserved higher-order physical residuals supply
supplementary evidence on real fields. Absolute error bounds for the complete
original-source continuum problem remain unavailable.

The prior independently evaluated production orders remain unchanged, with their
source hash recorded in existing-physical-gate.json. These are **not** presented
as a fresh physical audit:

| Band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| inner | H | 1.740849102 | 1.929826071 | 1.996174863 |
| inner | M2 | 0.330580315 | 0.855566995 | 1.112652071 |
| outer | H | 1.745667434 | 1.963989971 | 2.018377225 |
| outer | M2 | 1.852947017 | 1.970457203 | 1.992557652 |

Result remains **9/12 FAIL**. No band, row, order or threshold is removed.

## Decision and reproduction

The evidence supports the existing direct-Hessian closure as an accuracy
improvement, already executed at all four production resolutions. It gives no
new accuracy-based reason to modify that closure so its fields counterbalance
the audit defect. Adding that defect to the solver, adjusting physical forcing,
weakening the test, or selecting coefficients from the observed ratios would
not be a justified accuracy correction here. Accordingly no new production
candidate is invented and no redundant production executions are claimed.

This is the scientific blocker permitted by the task: no defensible **next**
correction satisfying the retained finite-grid gate is established. The gate's
interpretation as a universal accuracy validator is contradicted by exact-field
counterexamples. This does not prove that a consistent 12/12 method cannot exist.

All **44 tests pass**, including the new independent 50-digit forcing test.
All **505 pre-existing files, including 153 NPZ archives, are hash-identical**.
Original source files, archives, checkpoints, solvers and frozen audit are intact.
New evidence is under `results/first-row-consistency-20261103`: environment,
sparse coefficient certificates, all manufactured arrays, physical row residual
arrays, GMRES histories, logs and SHA256SUMS. The initial radial-energy import
failure and corrected execution logs are both preserved.

Replay from repository root with a fresh tag:

    DEE_ROW_TAG=unique-first-row-replay bash research/reproduce_run43_first_row.sh

Versions, including mpmath 1.3.0, are in environment.json. Output paths refuse
overwriting. The replay runs only these new checks, manufactured solves and
regression tests; it never changes or re-solves a production source problem.

No merge, promotion, acceptance change, physical DEE validation claim or
Physical Review D submission is made.
