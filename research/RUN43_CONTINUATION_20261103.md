# Run 43 continuation: original sources and radial closure

This study is isolated on `research/run43-original-continuation-20261103`,
branched from `7e8e6bac4bbcbf424465a7b22252744ec672ffa9`.
The remote branch retains the entire prior tree and ancestry. The local checkout
was recovered through GitHub blob reads after the configured execution proxy
refused connections. Its local initial commit is explicitly a partial recovered
source snapshot, not a clone or a replacement for the remote ancestry.
`results/continuation-20261103/recovery-provenance.json` records file blob hashes.
The November 2 successful source, settings, execution log and independent audit
are copied under `preserved-nov02/`; all original remote artifacts remain intact.

## Continuum equations and discrete distinction

The original physical sources are fixed Run38 fields and normal velocities.
Their energy is recomputed with the physical inverse metric:

```
E(psi) = temporal/2 + spatial * psi^(-4)/2 + U
S_i = -K_AB Pi_A D_i phi_B
Delta psi + |LW|^2 psi^(-7)/8 + 2 pi G E(psi) psi^5 = 0
partial_j (LW)_ij = 8 pi G psi^6 S_i
(LW)_ij = partial_i W_j + partial_j W_i - 2 delta_ij div(W)/3
```

Here `G=0.001`, `gamma_ij=psi^4 delta_ij`, maximal slicing `K=0`,
radial Dirichlet data `psi=1,W=0` at r=0.05,0.25, and periodic angles
are unchanged. The physical momentum residual is
`psi^(-10) partial_j(LW)_ij - 8 pi G psi^(-4) S_i`.

The original expanded discrete operator is
`E_i=L W_i + D_i sum_j D_j W_j/3`.
The experimental tensor operator is
`T_i=sum_j D_j(D_i W_j + D_j W_i - 2 delta_ij sum_k D_k W_k/3)`.
Subtracting gives the exact discrete identity

```
T_i-E_i = (sum_j D_j D_j - L) W_i + sum_j [D_j,D_i] W_j.
```

The test computes both terms independently on random vectors. In smooth
Cartesian calculus, derivatives commute and `sum partial_j^2=Delta`, so the
operators represent the same continuum constraint. They are not identical
finite-grid equations. The scalar conservative face-flux operator also has the
same smooth limit as the expanded chart Laplacian. Both scalar and vector
stabilizations subtract `kappa delta_angle^4/(h_angle^2 scale^2)` with fixed
`kappa=1/4`. Taylor expansion gives an O(h_angle^2) smooth-field correction;
its Fourier symbol is nonpositive. This is consistency and tested-mode damping,
not a uniform stability proof or a guarantee of convergence for arbitrary data.

`compare_run43_operators.py` verifies the identity, compares smooth manufactured
truncation errors and checks the original source against direct 4D stress-energy
projections. Manufactured fields are used only for operator testing.

## Justified numerical correction

The original radial endpoint derivative has leading error `-h^2 f'''/3`,
where the adjacent centered derivative has `+h^2 f'''/6`. At the first interior
row, composing these produces `D(Df)-f'' = h f'''/4 + O(h^2)` at the inner
boundary and the opposite sign at the outer boundary. Thus second-order
first derivatives need not give second-order composed second derivatives there.

The corrected four-point endpoint moment conditions are
`sum w=0, sum w*j=1, sum w*j^2=0, sum w*j^3=1`, with offsets measured
from each endpoint and physical weights divided by h. They match the centered
`+h^2 f'''/6` term. For smooth fields the mismatch then contributes O(h^2),
including when differentiated at the endpoint itself. Geometry multiplication
still gives ordinary truncation commutators; it is not assumed to commute
exactly at finite resolution.

`EndpointGrid` changes only the radial derivative used in longitudinal tensor
and tensor divergence. The scalar face-flux stencil, source, stabilization,
nonlinear and linear tolerances and frozen audit remain unchanged. It also
checks the true unpreconditioned vector residual against the existing GMRES
tolerance, rejecting nonfinite residuals. Its executed correction source is
preserved alongside the baseline solver snapshot.

This candidate is tested separately from the unchanged successful November 2
discretization. A better truncation closure need not pass a coarse-grid residual
ratio, and no failed check is removed or retuned.

## Regression repair and reproduction

All 32 recovered tests passed before changes; no current test failure was
fabricated. The stale regression workflow ran only one Hessian test, targeted
an older branch, and depended on historical field files. It now discovers the
complete suite, checks the frozen audit hash, and includes the current branch
and matrix-free files in its triggers. The tensor tests isolate experiment
environment variables; two new tests check radial composition convergence and
exact preservation of matter and scalar operators. Both ordinary and production
experiment environments are tested. Full before/after logs are retained.

Run from the repository root, with Python and NumPy/SciPy versions recorded in
the environment evidence. Choose fresh output names; existing evidence must
never be overwritten.

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python research/compare_run43_operators.py operator-comparison.json
OPENBLAS_NUM_THREADS=1 DEE_EXPERIMENT=new-baseline DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier python -u research/run43_consistent_momentum.py
OPENBLAS_NUM_THREADS=1 DEE_EXPERIMENT=new-endpoint DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier python -u research/run43_endpoint_continuation.py
OPENBLAS_NUM_THREADS=1 python research/audit_run43_arrays.py research/results/new-baseline
OPENBLAS_NUM_THREADS=1 python research/audit_run43_arrays.py research/results/new-endpoint
OPENBLAS_NUM_THREADS=1 python research/validate_original_physics.py research/results/new-baseline
OPENBLAS_NUM_THREADS=1 python research/validate_original_physics.py research/results/new-endpoint
```

The independent array auditor reloads and hashes fields and verifies the frozen
Git blob `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`. Its twelve orders cover
H and M2 in both fixed physical boundary bands for 8→16,16→32,32→64.
Every order must be finite and at least 1.5. The supplementary original-physics
audit uses different radial and angular derivatives, direct source projections,
and nested-grid solution differences; it does not replace these criteria.

These are finite-shell fixed-source initial constraints. Neither the original
nor experimental operator solves the full Einstein–DEE evolution and matter
field equations. No general nonaxisymmetric validation, physical DEE proof,
production promotion, merge, or journal submission is claimed or performed.

## Executed results
Both studies completed all original-source resolutions n=8,16,32,64. The separate-process array audit reports **12/12 PASS** for the unchanged November 2 discretization and **8/12 FAIL** for the matched endpoint candidate. The source and audit were not altered to obtain acceptance.

### Unchanged successful discretization
| Band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| inner | H | 1.740847585 | 1.929825463 | 1.996175084 |
| inner | M2 | 1.875425682 | 1.970278594 | 1.994731881 |
| outer | H | 1.745663903 | 1.963989315 | 2.018376693 |
| outer | M2 | 2.229142986 | 2.132462281 | 2.071256761 |

Nested-grid volume-weighted solution orders: {"psi": [1.9520561024297924, 1.9972615377554641], "W": [1.8973977020113797, 2.032459072833345]}. These are self-convergence estimates, not absolute solution-error bounds.

Supplementary inner M orders: [-0.16187043048169483, 0.6964093398200806, 1.1301952819022374]. This independent derivative diagnostic has no substituted acceptance threshold.

### Matched radial endpoint candidate (rejected)
| Band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| inner | H | 1.740847590 | 1.929825391 | 1.996175110 |
| inner | M2 | 0.234124077 | 0.727095846 | 1.116755252 |
| outer | H | 1.745663928 | 1.963989215 | 2.018376479 |
| outer | M2 | 1.998783011 | 1.532308823 | 1.394913925 |

Nested-grid volume-weighted solution orders: {"psi": [1.9520560070913506, 1.9972615993390466], "W": [1.903684177202203, 1.9985219676748534]}. These are self-convergence estimates, not absolute solution-error bounds.

Supplementary inner M orders: [0.40071163418974537, 1.2315868087591662, 1.764872449740053]. This independent derivative diagnostic has no substituted acceptance threshold.

### Scientific interpretation and preservation
The successful frozen gate is reproduced, including bitwise-identical compressed numerical field files at all four resolutions relative to the preserved November 2 Git blobs. Source-projection differences are zero for the analytic reconstruction; direct 4D sample discrepancies are at most 1.34e-15. All solutions are finite and positive with exact radial Dirichlet values. Full GMRES histories and every completed nonlinear checkpoint are retained.

The radial closure is analytically justified and its scalar composed derivative tests approach second order at endpoints and adjacent rows. Nevertheless, its frozen M2 inner-band orders are 0.2341,0.7271,1.1168 and its final outer M2 order is 1.3949. Those four scientific failures are preserved. The unchanged audit uses the original endpoint closure; its derivative mismatch remains visible on the corrected fields. Better smooth-field truncation does not imply this frozen coarse-grid gate will pass. No tuning, row exclusion or promotion follows.

The successful baseline also has supplementary inner M orders -0.1619,0.6964,1.1302 despite its 12/12 frozen pass. This is a genuine limit on claiming uniform boundary accuracy. The endpoint candidate improves those supplementary orders to 0.4007,1.2316,1.7649 but does not meet the frozen gate. A closure satisfying both robust boundary accuracy and the frozen twelve checks is not established by these experiments. No numerical breakdown prevented completion of either study.

Thirty-four tests pass in both default and production experiment environments. The replay shell script passes bash syntax validation. The local portable evidence archive and its SHA-256 accompany the complete source-hashed fields/checkpoints/logs; the separate remote research branch preserves the original ancestry and adds evidence without touching baseline branches.
