# Experimental consistent momentum study

Base: `6cc7364b181a6c2bc7daa7b9489d82359d09d087` on
`research/run43-inner-boundary-repair-20261101`. Work is isolated on
`research/run43-consistent-momentum-20261102`; no promotion is authorized by this report.

## Discrete discrepancy

Let D_i denote the collocated Cartesian first derivatives and L the direct chart
Laplacian. The original momentum solve uses L W_i + D_i D_j W_j / 3.
The independent tensor route uses T_i = D_j(D_i W_j + D_j W_i - 2 delta_ij D_k W_k/3).
Their difference is `(sum_j D_j D_j - L) W_i + sum_j [D_j,D_i] W_j`.
Thus replacing only the endpoint derivative cannot remove the entire discrepancy.
The existing three-point endpoint closure additionally makes the first interior
row of the composed radial second derivative only first-order accurate.

The experimental matrix-free solve evaluates T directly. Centered collocated
first derivatives underresolve Nyquist modes. The candidate subtracts
`kappa * (delta_theta^4/(h_theta^2*r^2) + delta_phi^4/(h_phi^2*q^2))`,
with kappa fixed at 1/4 before the first coarse run. Each fourth difference is
nonnegative in Fourier space and the correction has a negative symbol. For smooth
fields its magnitude is O(h_angle^2); it changes the discretization, not the matter
source or continuum equations. No radial fourth difference is added.

## Scalar investigation

The first candidate retains the original scalar operator and fails inner H on
8-to-16 (1.48875886). A separate fourth-order scalar chart operator fails that
check more strongly (0.45269922); higher formal order alone does not guarantee a
coarse-grid residual ratio. Both attempts are retained.

The next candidate uses conservative scalar face fluxes, matching the algebra of
the unchanged frozen Hamiltonian diagnostic, with the same O(h_angle^2) negative
fourth-difference stabilization. Its audit residual measures the stabilization
defect rather than the expanded-versus-flux mismatch. This is a stabilized
second-order discretization, not an exact enforcement or alteration of the audit.
The resulting gate is a constraint-residual convergence test; it cannot establish
solution error convergence by itself. No parameter was fitted to convergence data.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 DEE_EXPERIMENT=new-study DEE_SCALAR_ORDER=flux \
  DEE_PRECONDITIONER=fourier python research/run43_consistent_momentum.py
python research/audit_run43_arrays.py research/results/new-study
```

Use a fresh experiment directory for each execution. `fields_N.npz` contains psi
and Cartesian W at `(N+1,N,N)`. Per-resolution JSON contains SHA-256 hashes,
unpreconditioned linear residuals, complete GMRES histories and all audit regions.
Checkpoints retain every completed nonlinear iteration. The full flux study also
records its executed source and numerical-library settings. Early exploratory
runs preceded source-snapshot instrumentation; their source evolution is described
here, and they must not be confused with the fully recorded four-resolution run.

The scalar Fourier-block preconditioner exploits coefficient independence of phi,
but keeps all Fourier modes; it imposes no symmetry reduction on the unknowns.
Newton solves use the same Jacobian and positivity checks, with independently
checked unpreconditioned linear residuals. Vector GMRES retains rtol=1e-10,
atol=1e-13, restart=60 and at most 100 restart cycles.

## Scope and limitations

The frozen audit file must retain Git blob
`5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`. Exactly twelve checks cover H and M2,
both physical boundary bands, and every pair of 8,16,32,64. Every order must be
finite and at least 1.5. A separate process reloads and hashes every solution and
recomputes the checks using that file.

Even a 12/12 result is not a proof of second-order continuum solution accuracy:
the composed radial endpoint defect remains first order locally. The checkerboard
tests establish damping of a tested angular mode and absence of an exact null mode
on one small full grid, not a uniform stability theorem for every resolution.
The independent audit is independent evaluation code, not an independent research
group's replication. Axisymmetric matter on a finite shell remains the benchmark;
there is no general nonaxisymmetric evolution, exterior matching, physical DEE
proof, production promotion or journal submission.

## Completed frozen-gate results

All four solutions were computed fresh with the same executed solver snapshot,
settings and matter source. A separate array-audit process found **12/12 PASS**.
The original frozen file is unchanged. Orders are:

| Physical band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| Inner [0.05,0.10] | H | 1.740848 | 1.929825 | 1.996175 |
| Inner [0.05,0.10] | M2 | 1.875426 | 1.970279 | 1.994732 |
| Outer [0.20,0.25] | H | 1.745664 | 1.963989 | 2.018377 |
| Outer [0.20,0.25] | M2 | 2.229143 | 2.132462 | 2.071257 |

The minimum is 1.7408475852, above the unchanged 1.5 threshold. All psi values
are finite and positive (1 to 1.000167023), radial Dirichlet values are exact,
and maximum algebraic momentum residuals are below 1.6e-12. The legacy diagnostic
key `expanded_momentum_solver_max_abs` here refers to the replacement stabilized
operator, not the original expanded operator. The 64-grid solve took 327 seconds.

`results/stabilized-flux-full/independent-audit.json` gives every unrounded order,
RMS residual, raw-array SHA-256, and independent fourth-order supplementary audit.
`results/consistent-momentum/independent-audit.json` preserves the **11/12 FAIL**
from the scalar-baseline candidate, including its failed coarse pair. The
fourth-order scalar exploratory run preserves its **3/4 FAIL** partial study;
it was rejected before spending resources on its fine grids.

For any scalar candidate solved to algebraic tolerance, the frozen H residual is
`-8*psi^-5*(L_flux-L_candidate)*psi` plus the algebraic nonlinear residual.
Changing L_candidate changes this truncation defect. Formal asymptotic stencil
order does not bound the first coarse-grid ratio: its leading coefficient can
cancel other chart-direction terms, and each fixed band gains radial samples on
refinement. This explains why the fourth-order scalar candidate did not repair
the coarse frozen check. The conservative stabilized candidate removes the
uncontrolled expanded-versus-flux defect and leaves a specified O(h_angle^2)
term. Passing this criterion does not erase the radial T truncation limitation.

All 32 unit tests passed, including independent tensor-route equality, angular
checkerboard damping, a full small-grid singular-spectrum null-mode check,
Fourier inverse accuracy and independent scalar flux consistency.

GitHub Actions run [41386724276](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/actions/runs/41386724276)
completed successfully on commit `4e3d4d5`, including fresh 8/16/32/64 solves,
32 tests, separate raw-array twelve-check evaluation, and evidence upload.
Its full decoded job log is preserved in `results/actions-41386724276.log`.
The workflow artifact holds the independent runner's raw arrays and checkpoints.
A local attempt to download that artifact received HTTP 403 from the signed Azure
storage redirect; the successful job log was retrieved through the GitHub connector.
Local arrays and checkpoints are fully preserved in this branch independently.

No production promotion or claim of comprehensive physical validation follows
from the numerical gate. The remaining technical research question is uniform
continuum accuracy and radial checkerboard stability of the endpoint closure.
No user intervention is required to reproduce the present successful gate.
