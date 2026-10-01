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
