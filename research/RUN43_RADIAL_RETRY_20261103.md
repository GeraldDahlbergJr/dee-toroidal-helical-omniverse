# Run 43 radial accuracy retry

Work is isolated on `research/run43-radial-consistency-retry-20261103`, based on
remote commit `20c25a10aa34538a948be83b56783b8c3fdfe5ca`. The original **12/12
passing baseline remains unchanged and accepted by its frozen numerical gate**.
Neither candidate below replaces it. No source, physical boundary data, frozen
criterion, prior field, checkpoint or historical result was overwritten.

## New numerical work

The previous matched endpoint experiment changed the first-derivative closure
but left collocated composed radial derivatives in the vector principal part.
Two new candidates were defined analytically, without fitting coefficients to
the frozen orders, and each was executed at n=8,16,32,64 with original sources.

1. **Radial principal-part correction.** Replace the composed radial second
   derivative by a direct three-point second derivative using
   `(I+er er^T/3)(Drr-Dr Dr)W`. The coefficient is fixed by the continuum
   principal symbol. This removes the first-order first-interior-row defect and
   supplies the direct radial elliptic response to checkerboard modes. The
   scalar face-flux operator and its fixed angular correction remain unchanged.
2. **Fourth-order direct Cartesian Hessian vector operator.** Use five-point
   first and periodic second derivatives, six-point radial second-derivative
   closures, and analytic toroidal connection terms for the Cartesian Hessian.
   Form `Delta W + grad(div W)/3` from these direct Hessians, avoiding composed
   radial derivatives. This operator has no added vector stabilization; direct
   second derivatives resolve the Nyquist mode. The baseline stabilized
   conservative scalar operator is retained exactly.

Both use the original nonlinear source coupling, G, matter coefficients, fixed
normal velocities, inverse physical metric in E, radial Dirichlet values and
periodic angular data. The longitudinal tensor used by the fourth-order
Hamiltonian uses the fourth-order first derivatives; this changes numerical
differentiation, not the physical source or continuum equation.

The existing scalar Fourier-block preconditioner is retained. All Fourier modes
remain present. Linear tolerances are unchanged, and actual unpreconditioned
vector residuals are checked. Every completed nonlinear iteration is checkpointed.

## Accuracy and regression evidence

Manufactured tests are restricted to operator accuracy and never substituted
for original-source gate fields. The radial correction approaches second-order
accuracy at both first interior rows. The fourth-order candidate approaches
fourth order at endpoints and adjacent rows, and on an angularly varying
Cartesian polynomial vector with independently known exact vector Laplacian.
Small-grid complete singular spectra detect no exact null mode in tested scope;
this is not a proof of uniform stability on arbitrary grids.

All **41 tests pass** in both default and production experiment environments.
The fourth-order radial operator's finest-pair orders are 4.075 at the first
inner row, 3.902 at the first outer row, and 3.901 for the largest endpoint
error. The angular Cartesian polynomial order is 3.984.

On original-source solved fields at n=64, the supplementary independently
implemented seven-point-radial/Fourier-angular momentum diagnostic gives:

| Discretization | Inner M RMS | Outer M RMS | Inner M order 32→64 |
|---|---:|---:|---:|
| Preserved 12/12 baseline | 1.0657272e-3 | 4.7117458e-5 | 1.130195 |
| New radial principal correction | 9.8956616e-5 | 1.6097389e-5 | 1.768275 |
| New fourth-order vector | 5.3014787e-7 | 1.6379507e-9 | 2.908646 |

Thus the new fourth-order candidate reduces this inner supplementary residual
by about 2010 times at n=64. This is evidence under a different diagnostic,
not a reclassification of the baseline's frozen pass or independent replication
by another research group. The fourth-vector nested-grid W orders are 3.5723
and 4.7072; scalar psi orders remain 1.9521 and 1.9973. These are self-convergence
estimates, not absolute error bounds. Supplementary inner H orders remain
0.6861,1.3789,1.7553 because the second-order scalar operator is unchanged.
Uniform full-system accuracy is therefore not claimed.

## Frozen criterion and measured obstruction

Both new candidates pass **9/12**, failing the three inner M2 checks. Their
failures are preserved. The frozen audit blob remains
`5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`, with every required order finite
and at least 1.5 for acceptance. The matched-closure and higher-order diagnostics
are supplementary only; they do not replace any of those checks.

The original endpoint first derivative has `-h^2 f'''/3` error while the centered
interior derivative has `+h^2 f'''/6`. Composing them introduces an O(h) term
in the first interior second derivative. For a smooth solution, one row with
O(h) error contributes O(h^(3/2)) to a fixed-band RMS. Other angular truncation
terms can cancel or reinforce it on finite grids. This does not imply every
coarse pair attains an order of 1.5.

At n=64 the new fourth-vector frozen inner M2 RMS is 1.9997339e-4 while the
supplementary high-order RMS is 5.3014787e-7. The norm of their difference is
1.9954500e-4. About 50.6% of the frozen squared inner-band residual lies on the
first interior row. The radial-principal candidate has about 91.4% there.
These measurements isolate a large diagnostic/discretization discrepancy;
they do not prove that satisfying both accuracy and the frozen coarse-grid
gate is mathematically impossible.

**No full replacement satisfying both requirements has been established.**
The retry delivers a substantially more accurate tested vector candidate, but
it is not promoted because it fails three frozen checks. The user's original
12/12 baseline is preserved. No threshold, boundary band, row selection or
physical source was changed to manufacture acceptance.

## Evidence and reproduction

`results/radial-retry-20261103/` contains full execution/regression logs,
manufactured-accuracy reports, residual decompositions, hashes and provenance.
Each candidate directory contains all four numerical fields, per-resolution
reports, complete GMRES histories, nonlinear checkpoints and executed sources.
The fourth-vector repeat directory checks reproducibility with explicit
effective scalar/vector settings. In its inherited `settings.json`, .25 refers
to the baseline scalar correction; `effective-settings.json` correctly records
zero stabilization for the direct vector operator.

From the repository root, select unique output names to preserve existing data:

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier DEE_EXPERIMENT=new-principal python -u research/run43_radial_consistency.py
OPENBLAS_NUM_THREADS=1 DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier DEE_EXPERIMENT=new-fourth-vector python -u research/run43_fourth_vector.py
OPENBLAS_NUM_THREADS=1 python research/audit_run43_arrays.py research/results/new-principal
OPENBLAS_NUM_THREADS=1 python research/audit_run43_arrays.py research/results/new-fourth-vector
OPENBLAS_NUM_THREADS=1 python research/validate_original_physics.py research/results/new-principal
OPENBLAS_NUM_THREADS=1 python research/validate_original_physics.py research/results/new-fourth-vector
```

The frozen audit commands return exit status 1 for these recorded candidates
because of genuine scientific gate failures; all fields and logs remain
available. No software or Krylov failure prevented completion. These are
fixed-source conformal initial constraints on a finite shell, not full
Einstein–DEE evolution or a physical DEE validation. No production promotion,
merge or journal submission was performed.

## All frozen convergence orders (independent array reload)

### principal
| Band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| inner | H | 1.740847579 | 1.929825517 | 1.996175305 |
| inner | M2 | 0.177153787 | 0.677075645 | 1.080077051 |
| outer | H | 1.745663925 | 1.963989145 | 2.018376277 |
| outer | M2 | 1.757659753 | 1.676579446 | 1.587360219 |

Recorded 9/12; no failed order omitted.

### fourth-vector
| Band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| inner | H | 1.740849102 | 1.929826071 | 1.996174863 |
| inner | M2 | 0.330580315 | 0.855566995 | 1.112652071 |
| outer | H | 1.745667434 | 1.963989971 | 2.018377225 |
| outer | M2 | 1.852947017 | 1.970457203 | 1.992557652 |

Recorded 9/12; no failed order omitted.

### fourth-vector-repeat
| Band | Constraint | 8→16 | 16→32 | 32→64 |
|---|---|---:|---:|---:|
| inner | H | 1.740849102 | 1.929826071 | 1.996174863 |
| inner | M2 | 0.330580315 | 0.855566995 | 1.112652071 |
| outer | H | 1.745667434 | 1.963989971 | 2.018377225 |
| outer | M2 | 1.852947017 | 1.970457203 | 1.992557652 |

Recorded 9/12; no failed order omitted.

The repeated fourth-vector solve produces byte-identical numerical field archives at all four resolutions. Its separate-process frozen and supplementary audits reproduce the first execution.
