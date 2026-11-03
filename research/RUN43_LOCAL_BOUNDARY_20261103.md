# Run 43 localized radial boundary correction

Research branch: `research/run43-boundary-continuation-20261103`.
Parent: `12b2b2a27e09a38a304afb9dff521c814a219901`, which retains the
November 2 successful baseline and all subsequent investigations.
Frozen audit Git blob: `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.

## Equations and derivation

The benchmark is the conformally flat, maximal initial-data reduction of the
coupled Einstein–DEE constraints, with original fixed matter fields and normal
velocities. It is not the full dynamical Einstein–DEE system. In Cartesian
components, let A_ij = partial_i W_j + partial_j W_i -
2 delta_ij div(W)/3. The equations solved remain

    Delta psi + A_ij A_ij psi^-7/8 + 2 pi G E(psi) psi^5 = 0,
    partial_j A_ij = Delta W_i + partial_i div(W)/3
                  = 8 pi G psi^6 S_i,
    E(psi) = temporal/2 + spatial psi^-4/2 + U.

The shell has r in [.05,.25], periodic angles, psi=1 and W=0 at both
radial boundary planes. S, temporal, spatial, U, G and the matter coefficients
are inherited unchanged. No third-normal-derivative boundary condition follows
from these Dirichlet data or from the prescribed matter.

For the frozen centered interior D_r and three-point one-sided endpoint D_r,
direct substitution at first interior row 1 gives

    (D_r D_r f)_1 = (3f0 - 5f1 + f2 + f3)/(4h^2),
    (D_rr f)_1   = ( f0 - 2f1 + f2)/h^2,
    (D_r D_r - D_rr)f_1 = (-f0 + 3f1 - 3f2 + f3)/(4h^2)
                               = h f'''(r1)/4 + O(h^2).

At the last interior row the leading odd derivative changes sign. Because
e_r is independent of r, the radial second-derivative principal symbol of
Delta W + grad div W/3 is exactly I + e_r e_r^T/3. The new correction is

    C W = P_edge (I + e_r e_r^T/3)(D_rr-D_r D_r) W,

where P_edge is one only on the first and last interior radial planes.
This cancels the first-order radial composition defect without modifying any
other row. Direct D_rr is centered and second order at both selected planes;
it uses only known Dirichlet values and existing unknowns. C is O(h) on those
two planes for smooth fields, and tends to zero. There is no adjustable
coefficient. This is a local consistency repair, not a claim of fourth-order
accuracy, a full-system stability proof, or a means of enforcing the audit.
The prior global radial repair additionally changed every bulk radial row;
this experiment isolates the adjacent-plane contribution.

The candidate starts from the November 2 tensor discretization, retaining
the original fixed .25 angular stabilization and stabilized conservative
scalar face flux. It does not modify the already accurate fourth-order
candidate. All physical solutions use the original sources; analytic cubic
data occur only in a unit test of the derived stencil coefficient.

## Solver versus frozen audit

Writing G_i for chart-assembled Cartesian first derivatives and L for the
direct chart Laplacian, the original expanded solver E and frozen tensor T obey

    E_i W = L W_i + G_i sum_j G_j W_j/3,
    (T-E)_i W = (sum_j G_j G_j-L)W_i + sum_j [G_j,G_i]W_j.

Thus a radial correction cannot remove every angular, product-rule or mixed
finite-grid discrepancy. The new solver is B = T - Q_angle + C, where
Q_angle is the preserved angular stabilization. Its exact frozen residual is

    M_frozen = psi^-10(B W - b) + psi^-10(T W - B W),
    b = 8 pi G psi^6 S.

An algebraically converged solve can consequently fail the frozen convergence
gate. The fourth-order Hessian candidate already removes the radial defect;
the parent investigation measures its first-row algebraic residual near
1e-13 while its frozen first-row residual is near 1e-3. That evidence is retained
as historical evidence and is not presented as a new fourth-order execution.

The physical-band volume RMS of a single O(h) plane is asymptotically
O(h^1.5), but the leading coefficient need not be constant across n=8,16,32,64.
It is unjustified to infer a finite-pair order of at least 1.5 from that
asymptotic statement. No rows, bands or thresholds are changed here.

## Reproduction and evidence

Use a fresh DEE_EXPERIMENT directory; the inherited driver refuses overwrites.

    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
      DEE_EXPERIMENT=unique-local-boundary DEE_SCALAR_ORDER=flux \
      DEE_PRECONDITIONER=fourier python -u research/run43_local_boundary.py
    OPENBLAS_NUM_THREADS=1 python research/audit_run43_arrays.py \
      research/results/unique-local-boundary
    OPENBLAS_NUM_THREADS=1 DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier \
      python research/check_run43_local_identity.py \
      research/results/unique-local-boundary research/results/unique-local-identity
    OPENBLAS_NUM_THREADS=1 python research/validate_original_physics.py \
      research/results/unique-local-boundary

The frozen audit exits 1 for a failed gate; retain its output and continue
collecting supplementary diagnostics. These diagnostics never replace the
twelve acceptance checks. Raw fields, nonlinear checkpoints, executed solver
snapshots, GMRES histories, true residuals, settings and all per-grid audit
regions are retained under `results/boundary-local-20261103`.

No merge, baseline replacement, candidate promotion or journal submission is
authorized by this research experiment.

## Actual execution: rejected, 8/12

Fresh original-source coupled solves completed at n=8,16,32,64, each with two
nonlinear iterations. Separate-process raw-array re-evaluation using the frozen
audit returned **8/12 PASS, 4/12 FAIL** (exit 1 as expected for rejection).
The result is not promoted. No coefficient was changed after viewing results.

| Fixed physical band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| Inner [.05,.10] | H | 1.740847590 | 1.929825391 | 1.996175110 |
| Inner [.05,.10] | M2 | 0.234124077 FAIL | 0.727095846 FAIL | 1.116755252 FAIL |
| Outer [.20,.25] | H | 1.745663928 | 1.963989215 | 2.018376479 |
| Outer [.20,.25] | M2 | 1.998783011 | 1.532308823 | 1.394913922 FAIL |

All twelve checks use the unchanged finite-order threshold >=1.5. No band rows
were omitted. Volume-weighted RMS residuals from the same frozen audit:

| Band/constraint | n=8 | n=16 | n=32 | n=64 |
|---|---:|---:|---:|---:|
| Inner H | 3.587365160405e-04 | 1.073317582015e-04 | 2.817039356200e-05 | 7.061294583437e-06 |
| Inner M2 | 9.128115113783e-04 | 7.760732932423e-04 | 4.688404553381e-04 | 2.161962966908e-04 |
| Outer H | 1.687614137141e-04 | 5.032409227341e-05 | 1.289900720128e-05 | 3.183936683199e-06 |
| Outer M2 | 5.406323161317e-05 | 1.352721400373e-05 | 4.676677716939e-06 | 1.778387964836e-06 |

| n | Solve seconds | Largest vector true residual norm |
|---:|---:|---:|
| 8 | 0.160 | 9.176762788663e-12 |
| 16 | 0.854 | 3.751504390497e-11 |
| 32 | 13.993 | 1.072701879175e-10 |
| 64 | 232.721 | 3.462344652087e-10 |

Every GMRES info is zero and every explicitly recomputed true residual meets
max(1e-13,1e-10 ||b||). These are norms, not per-point residual maxima.
Checkpoints preserve both coupled iterations. Raw-array validation confirms
finite positive psi and exact original radial Dirichlet values.

The first-row residual decomposition is:

| n | Physical algebraic RMS | Operator mismatch RMS | Frozen M2 RMS |
|---:|---:|---:|---:|
| 8 | 5.729701128412e-14 | 1.350907869110e-03 | 1.350907869100e-03 |
| 16 | 2.616758319034e-14 | 1.760486808134e-03 | 1.760486808135e-03 |
| 32 | 1.172008254683e-14 | 1.561418100465e-03 | 1.561418100464e-03 |
| 64 | 2.746341481636e-15 | 1.038103330866e-03 | 1.038103330866e-03 |

The identity closes over all interior rows to 5.02e-18. Thus unresolved linear
solve error does not explain these failures. Restoring radial consistency changes
the solution and leaves the frozen first-order composition defect observable;
keeping the baseline angular stabilization does not compensate for it. This is
a measured residual explanation, not an absolute solution-error bound.

Supplementary original-action reconstruction matches energy and momentum sources
exactly, with 4D stress-tensor sample discrepancies <=1.34e-15. Its separate
higher-order diagnostic is retained, including endpoint values; no supplementary
order is substituted for acceptance. It does not establish full dynamical
Einstein–DEE validation or uniform stability of this discretization.

Fresh frozen re-evaluations of the *preserved* baseline and fourth-order raw
arrays reproduce 12/12 and 9/12 respectively. Their original files were not
rewritten. The full first-row sparse-matrix certificate also reproduces the
radial/angular/Leibniz identity at relative coefficient errors 5.17e-16 and
5.94e-16 for n=8,16. The new cubic test verifies the signed boundary correction
and unchanged bulk rows, matter source and scalar operator. All 45 regression
tests pass.

**Scientific blocker:** neither the fourth-order candidate nor this independently
justified local boundary repair satisfies the retained finite-grid gate. There
is no evidence here for a further principled correction that will achieve 12/12.
Consistency alone does not guarantee the requested finite-pair residual orders;
changing the physical sources, fitting a compensating audit defect, or imposing
an unsupported third-normal-derivative boundary condition is not justified.
This is not a proof that no consistent passing discretization can exist.
The preserved November 2 baseline remains the successful reference.

All 573 pre-existing files are verified byte-identical after execution. The
frozen audit blob remains 5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22. Evidence,
preservation manifests, environment, logs, source bundle and checksums are in
`results/boundary-local-20261103-evidence`. No merge, promotion or Physical
Review D submission was performed.
