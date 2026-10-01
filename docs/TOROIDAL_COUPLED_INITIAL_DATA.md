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
