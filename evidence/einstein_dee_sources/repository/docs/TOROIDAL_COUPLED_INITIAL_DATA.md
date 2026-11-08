# Coupled toroidal-shell initial-data benchmark

This extends frozen Run #38 (`2ff3b7ecbd01ad8e2aede2823864d6627275ab43`) by solving a conformal Hamiltonian equation and all three vector momentum equations on a finite toroidal shell. It preserves the benchmark matter fields; it does not freeze a physical exterior boundary model.

## Explicit construction

The flat background has line element `dr² + r² dtheta² + (R+r cos(theta))² dphi²`, with `R=1`, `r in [0.05,0.25]`, and periodic angles. Set `gamma_ij=psi^4 bar_gamma_ij`, `K=0`, and `A^ij=psi^-10 (bar L W)^ij`. The vector and tensor are represented in Cartesian components on the toroidal grid, so the Cartesian derivatives include the reciprocal toroidal basis.

Hold the Run #38 field values, spatial gradients, and normal velocities fixed. In particular, recompute `E = temporal/2 + psi^-4 spatial_background/2 + U` as the geometry changes, while the physical momentum covector is fixed. The mixed Theta–Psi kinetic term remains present.

Solve

- `bar_Delta psi + |bar L W|² psi^-7/8 + 2 pi G E psi^5 = 0`;
- `bar_Delta W^i + (1/3) bar_D^i bar_D_j W^j = 8 pi G psi^6 bar_gamma^ij S_j`.

Use `G=0.001`, `psi=1` and `W=0` at both radial boundaries. These are benchmark Dirichlet conditions on a shell, not asymptotic-flatness conditions or exterior matching. Newton iteration solves the scalar equation; preconditioned GMRES solves the vector equations; an outer iteration couples them.

## Validation and scope

Run `OPENBLAS_NUM_THREADS=1 python reproduce/solve_toroidal_helical_coupled_constraints.py` from the repository root. Refinement doubles all three grid spacings: `9×8×8`, `17×16×16`, `33×32×32`. The program saves JSON diagnostics and NPZ arrays containing the conformal factor and all three Cartesian vector-potential components.

Physical Hamiltonian residuals are reconstructed from scalar curvature, `K_ij K^ij`, and the metric-dependent matter energy. Momentum is independently reconstructed by taking the divergence of the longitudinal tensor, rather than reusing the expanded vector-equation residual. Their discrete stencils differ by truncation error. The independent momentum convergence gate uses the same fixed radial bulk interval `[0.10,0.20]` at every resolution; full-interior maximum residuals are also reported. A bulk gate does not establish boundary convergence or uniform constraint satisfaction on the whole shell.

The Run #38 stress is axisymmetric even though the phase is written `theta-3 phi`. Cartesian vector components have phi dependence through the basis. A three-dimensional grid does not establish general nonaxisymmetric validation. A subsequent benchmark must add genuine phi-dependent stress, check boundary-adjacent residual convergence, and check the initial data through a separate geometry implementation before making stronger claims.

Analytic manufactured Cartesian polynomial tests verify the vector differential operator under angular and radial refinement. Separate matter tests compare the source with the existing physical-metric projection implementation.

## Correction to the Hamiltonian-only precursor

The previous scalar stencil composed centered first derivatives and set radial derivative values to zero at endpoints. This permitted alternating grid values and produced a conformal-factor minimum near 0.0096 even for weak gravity. It is replaced with the ordinary nearest-neighbor second derivative plus the analytic toroidal first-derivative coefficients. The precursor now also recomputes matter energy with the physical inverse metric. Its momentum remains explicitly unsolved; it is retained as a precursor only.

## Measured local result

The three-grid bulk benchmark passed. Independent vector momentum RMS was 8.73680837e-4, 2.50880457e-4, and 6.57976734e-5, with observed orders 1.8001 and 1.9309. On the finest grid the physical Hamiltonian maximum was 1.3842e-10, and the expanded vector solver maximum was 2.8476e-13. Independent full-interior component maxima remained about 9e-4 near the radial boundaries. All 21 unittest-discovered tests passed, including the two new manufactured/source checks. These values support the scoped bulk benchmark; they do not justify a uniform whole-shell or general nonaxisymmetric claim.

## Frozen Run42 boundary audit

Run `python reproduce/diagnose_toroidal_boundary_convergence.py` to audit the
hash-verified Run42 archive without changing equations, matter seed, boundary
values, or solved fields. The audit reproduces the archived bulk norm and uses
fixed physical radial bands. Its primary gate requires both independent
face-flux Hamiltonian and second-order tensor-divergence momentum RMS to have
observed order at least 1.5 on both refinement pairs in each boundary band.

The result is **BOUNDARY_GATE_NOT_MET**. The inner band `[0.05,0.10]` has
Hamiltonian orders 1.4888 and 1.7354, and momentum orders 0.4528 and 0.8215.
The outer band `[0.20,0.25]` passes the declared criterion: Hamiltonian orders
2.1439 and 2.0575, momentum orders 1.8450 and 1.9681. Bulk momentum retains
orders 1.8001 and 1.9309. Radial endpoints are reported separately because
Dirichlet values alone do not enforce endpoint PDE compatibility.

A supplementary fourth-order derivative audit gives inner-band momentum orders
2.7636 and 2.3223. This suggests derivative truncation contributes to the inner
boundary behavior; it does not correct the second-order solved fields or replace
the primary gate. The next investigation should isolate the inner boundary
operator with manufactured solutions before changing the solver. A green CI
run means the diagnostic completed and passed integrity checks; the JSON's
scientific status remains explicitly unresolved.

## Inner operator investigation during Run44

`python reproduce/investigate_toroidal_inner_operator.py` records manufactured
fields, six radial refinements, and hash-verified Run42 field comparisons in
`toroidal_inner_operator_investigation.json`. No production solver is changed.

The physical inner-band momentum residual is reproduced by the difference
between the independently differentiated longitudinal tensor and the expanded
vector operator: RMS 1.0160549e-3, 7.4236816e-4, and 4.2008194e-4. This points to
operator consistency, rather than insufficient convergence of the algebraic
vector solve. The difference decomposes into a composed-versus-direct scalar
Laplacian defect and noncommuting discrete Cartesian derivatives. Their vector
norms must not be added as scalars; cancellation is present.

At the first interior radial row, composing centered first derivatives uses an
endpoint derivative whose leading truncation term differs from the interior.
The usual endpoint formula has error `-h² f'''/3`, versus `+h² f'''/6` for the
centered formula. Differentiating that error jump produces an O(h) term. Even a
third-order endpoint derivative leaves a mismatch with the centered interior.
A four-point audit closure matching the centered leading term restores roughly
second-order composition in the isolated radial test: the last refinement pair
has order 2.0193, versus 1.0738 for the original closure. In the full manufactured
vector test the matched closure gives first-inner-row orders 1.8400, 1.9130,
and 1.9533.

This is not yet a remedy for the frozen solution. Changing only the endpoint
formula in its audit improves inner-band orders from 0.4528/0.8215 to
0.8136/1.3303, still below the declared 1.5 gate. The production expanded
operator itself also uses composed derivatives in `grad(div W)`. A consistent
candidate vector operator needs a new solve and independent boundary study;
relabeling the existing result is not justified. Endpoint compatibility,
angular resolution near the small tube radius, and stability remain to be
checked. Two additional tests verify cubic composition and reproduce the
production vector operator independently; all 25 local unit tests pass.

## Fresh matched-endpoint candidate: rejected as a standalone remedy

`OPENBLAS_NUM_THREADS=1 python reproduce/solve_toroidal_matched_endpoint_candidate.py`
performs fresh solves on `9x8x8`, `17x16x16`, and `33x32x32`. A separate subclass
changes only both radial first-derivative endpoint closures and reconstructs
the vector operator. The original production solver, equations, seed, geometry,
boundary values, tolerances, and independent acceptance criterion are preserved.
Candidate NPZ fields and their hashes are recorded separately.

The result is **CANDIDATE_BOUNDARY_GATE_NOT_MET**. The original independent
inner-band momentum orders are 0.452768 and 0.821458, essentially unchanged.
The supplementary matched-endpoint audit gives 0.813589 and 1.330309, also below
the 1.5 threshold. Inner-band independent Hamiltonian orders remain 1.488759
and 1.735411; the outer band passes. Conformal-factor changes are at roundoff,
and the maximum Cartesian vector-potential change is below 4.7e-10 across all
three grids. The algebraic solves converged; this candidate fails the scientific
boundary gate and is not promoted, frozen, or tagged. All 26 local tests pass,
including unchanged matter/scalar-interior checks and an independent comparison
of the candidate expanded operator.

The isolated endpoint correction is therefore insufficient for these fields.
The next investigation should test a vector discretization with compatible
second derivatives, distinguish scalar Laplacian and discrete-commutator
contributions, and separate radial from angular refinement. Any proposed remedy
must be checked through a fresh solve and an independent audit rather than
changing the acceptance criterion or merely reducing the solver tolerance.

## Direct Hessian candidate and directional derivative audit

`solve_toroidal_direct_hessian_candidate.py` forms `grad(div W)` from direct
chart second derivatives, mixed derivatives, and the toroidal connection,
then rotates the Hessian to Cartesian components. It retains the scalar solver,
Run42 continuum equations, matter seed, geometry, boundaries and tolerances.
The production solver remains unchanged. Manufactured inner-first-row and
inner-band operator orders are 2.2755 and 2.1099 on the 9/17 radial pair.
An independent Cartesian polynomial test also converges; all 28 tests pass.

Three fresh candidate solves completed, but the original boundary criterion
still fails: inner-band momentum orders 0.50824 and 0.85244, versus Run42's
0.45277 and 0.82146. Independent flux-Hamiltonian orders remain 1.48876 and
1.73541. Outer-band criteria pass. Fourth-order momentum auditing is
supplementary: inner orders 2.62185 and 3.06471 do not replace the primary gate.
The candidate fields, hashes and diagnostics are saved separately and are not
promoted, frozen, or tagged.

`diagnose_toroidal_derivative_directions.py` isolates differentiation order in
the saved fields without changing or interpolating them. For frozen Run42,
using fourth-order radial derivatives with second-order angular derivatives
gives inner-band momentum orders 1.18688 and 1.78079. Fourth-order angular
derivatives with second-order radial derivatives give 0.39066 and 0.76111.
The direct-Hessian candidate shows the same pattern. This supports radial
differentiation as a major contributor to the measured inner-boundary problem;
it is not independent radial/angular *grid refinement* or proof of a remedy.

Before another solver change, the next useful check is controlled radial
refinement with angular resolution fixed, followed by the complementary angular
refinement, using an independent boundary derivative implementation. The
Hamiltonian coarse-pair limitation also needs investigation. Lowering solver
tolerances or substituting the supplementary audit would not establish the
original boundary gate.
