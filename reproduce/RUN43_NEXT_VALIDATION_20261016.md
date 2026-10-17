# Run43 next validation decision — 2026-10-16

## Evidence now verified

The independent manufactured polynomial test completed in GitHub Actions run 39000842272, with second-order inner-band observed orders 1.324, 1.810, 1.942 and fourth-order 3.141, 3.871, 3.975. This demonstrates asymptotic accuracy for that particular smooth polynomial field; it does **not** establish correctness of the physical solution.

The frozen physical 33→65 inner-boundary momentum order remains 1.092 (< 1.5 gate). The solver's low algebraic residual is not a substitute for this independent physical constraint.

## Next controlled experiment

1. Recompute the *same* independent momentum residual on archived 9/17/33/65 candidate fields, preserving original SHA256 hashes and unmodified original metrics.
2. Produce radial-shell residual profiles for each resolution (especially first and second interior radial rows). Record per-shell weighted RMS and maximum; distinguish first-row contamination from error distributed across the fixed inner band.
3. Compare the numerical momentum operator used in the solver (direct Hessian) with the independent double-gradient divergence at the same saved fields, and test the difference on a smooth manufactured field that satisfies the vector boundary values.
4. Audit whether the vector boundary values and matter source imply a compatible endpoint momentum equation. Exact Dirichlet values alone do not prove PDE compatibility.
5. Only after diagnosing the discrepancy, propose a change on this branch and rerun coarse grids before a new 65-grid solve.

Do not relax the frozen >=1.5 convergence criterion, alter original checkpoints, or present the manufactured 1.942 as Run43 physical convergence. This note is an investigation plan, not a new solver result.
