# Run 43 streamed manufactured qualification

**Overall manufactured qualification: NOT YET ESTABLISHED.** Original-source validation and independent physical reproduction were gated. PRD readiness remains on hold.

## Memory repair and equivalence

The pinned run failed in the generated symbolic momentum forcing, at the x88 expression requesting another 54 MiB array on a 192³ quadrature grid. Sparse coordinate inputs broadcast to whole-domain common-subexpression temporaries while analytic values and gradients remained live. The last recorded historical RSS pertains to n=16. The failed n=32 process has no recorded peak RSS.

The repair streams eight radial quadrature planes at a time. It retains the four forcing components, evaluates the same symbolic expressions and quadrature nodes, and accumulates the same weak forms. Angular transpose contractions precede the radial transpose so a slab is not expanded before angular contraction. Scalar and vector applications, loads and strain-dependent terms use the same slabs. The equations, physical constants, spline space, quadrature order, and boundary prescriptions are unchanged.

At n=8 and n=16, forcing was bitwise identical. The largest relative difference in seven tested forcing/operator/load comparisons was 2.62e-16 under rtol 5e-12 and atol 5e-14. The two completed historical field and residual metric sets agree within 2.87e-15 absolute, with exact trace agreement. The n=32 preflight reached 532620 KiB peak RSS, with forcing taking 10.19 s and one vector application 4.19 s. See equivalence.json, results_equivalence.json and profile_memory.jsonl.

## Authorized rerun

| n | Cumulative seconds | Peak RSS MiB | Max linear relative residual | Final scalar / vector relative defects | Exact q6 Schur upper bound | Largest quadrature proxy / RMS |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 13.50 | 111.18 | 1.50399e-12 | 3.43025e-16 / 1.25128e-12 | 0.00107877 | 2.51463e-05 |
| 16 | 111.92 | 174.52 | 1.87626e-12 | 6.35558e-16 / 1.87593e-12 | 0.000817673 | 1.51076e-06 |
| 32 | 848.25 | 640.79 | 1.52492e-12 | 1.46024e-15 / 1.50534e-12 | 0.000781947 | 1.26132e-07 |

Driver process exit code: 0. Driver diagnostic status: INCONCLUSIVE. Stop: all grids complete; final uncertainty and quadrature review required.

| Grid pair | Min order | Max order | Region-component checks |
|---|---:|---:|---:|
| 8→16 | 2.963555124 | 5.970288311 | 20 |
| 16→32 | 2.953881763 | 4.415878813 | 20 |

Orders use the frozen q6/q8 quadrature sensitivity proxy. The 20 metrics per pair cover H and all three Cartesian momentum components in both boundary bands, bulk and both first interior planes. See results.json for the full rates and residuals.

| n | Analytic field RMS errors (psi, Wx, Wy, Wz) |
|---:|---|
| 8 | 2.912059e-07, 4.094535e-07, 2.598128e-06, 4.386608e-06 |
| 16 | 4.22189e-09, 4.585746e-09, 2.856981e-08, 4.308036e-08 |
| 32 | 1.049712e-10, 1.082961e-10, 6.385593e-10, 9.21209e-10 |

## Acceptance matrix

| Criterion | Verdict | Scope |
|---|---|---|
| Forcing and weak-kernel equivalence | PASS | Bitwise forcing and roundoff-scale weak-form differences on both overlap grids. |
| Complete one-thread, 4 GiB, 30-minute rerun | PASS | 3/3 completed grids. |
| Measured strong-residual rates ≥1.5 with quadrature proxy | PASS | 40/40 measured orders. This row is limited to the recorded proxy. |
| Relative linear residual ≤1e-10 | PASS | Explicit residual recomputation. |
| Radial Dirichlet and angular periodic traces | PASS | Zero radial samples. Periodicity follows from translated cardinal splines with modulo coefficient indexing. The seam samples also report zero, but wrap to the same argument. |
| No unexpected coupled discrete nullspace | PASS | Exact Gauss-6 Jacobian at archived coefficients. See CONDITIONING.md. |
| Analytic field-error comparison | PASS | Independent piecewise-polynomial field evaluation against the analytic fixture. RMS and sampled maxima archived. |
| Justified total evaluator uncertainty ≤10% | NOT YET ESTABLISHED | Quadrature sensitivity excludes a certified truncation allowance and real-node evaluation error. Final algebraic defect and inverse bounds alone do not enclose propagated assembly and nonlinear solver effects. |
| Original-source constraints and independent reproduction | NOT YET ESTABLISHED | Qualification prerequisite did not pass. No original-source solve or n=64 run. |

## Exact remaining evidence

For each archived field, component and region, enclose the weighted squared-residual integral and volume normalization, including rounding at quadrature nodes and quadrature truncation. One concrete test is outward-rounded cellwise analytic bounds for each coordinate with tensor Gauss error propagation, including angular-only quadrature on the first interior planes. Relative enclosures of 15% for the numerator and 3% for the denominator suffice for a 10% RMS allowance. These are sufficient bounds, not new acceptance thresholds. Where the protocol compares to the exact Galerkin root, bound coefficient perturbation from assembly quadrature and final nonlinear termination, then propagate that perturbation to the residual metric. The q6 Jacobian certificate can support that analysis but is not a nonlinear root enclosure.

No new DEE solve is required to audit the fixed archived fields. The authorized repair and single qualification attempt are spent. No further repair, original-source experiment, or manuscript preparation was initiated.
