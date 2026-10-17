# Run #43: preserved boundary convergence diagnostic

Read-only cross-resolution audit. No solver, fields, checkpoints, or production workflows changed. Source: archived candidate JSON at `coupled-einstein-dee-evolution` and Run #5396 artifact `11847550871` (commit `be21c4a`). Values are volume-weighted RMS in fixed physical boundary bands. Refinement ratio is 2; observed order = log2(error_coarse/error_fine).

| Region | Audit | 9×8×8 | 17×16×16 | 33×32×32 | 65×64×64 | p(9→17) | p(17→33) | p(33→65) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| inner boundary band | Hamiltonian flux | 1.0026e-3 | 3.5725e-4 | 1.0729e-4 | 2.9412e-5 | 1.489 | 1.735 | 1.867 |
| inner boundary band | Momentum order 2 | 1.0894e-3 | 7.6595e-4 | 4.2422e-4 | 1.9897e-4 | 0.508 | 0.852 | 1.092 |
| inner boundary band | Momentum order 4 | 8.8840e-4 | 1.4433e-4 | 1.7250e-5 | 7.0352e-6 | 2.622 | 3.065 | 1.294 |
| outer boundary band | Hamiltonian flux | 2.9208e-5 | 6.6090e-6 | 1.5877e-6 | 3.9016e-7 | 2.144 | 2.058 | 2.025 |
| outer boundary band | Momentum order 2 | 9.1317e-4 | 2.5429e-4 | 6.5020e-5 | 1.6350e-5 | 1.844 | 1.968 | 1.992 |
| outer boundary band | Momentum order 4 | 4.6204e-4 | 4.2117e-5 | 4.5406e-6 | 7.4451e-7 | 3.456 | 3.213 | 2.609 |

## Assessment

The Run42 boundary gate requires second-order Hamiltonian and momentum RMS observed orders ≥1.5 on both refinement pairs in each boundary band. The candidate's 65-grid refinement provides an additional pair, but the second-order inner-boundary momentum audit is below 1.5 on the finest pair. The fourth-order inner momentum audit also falls below 1.5 on that pair. Thus **boundary validation is NOT established**, even though the 65-grid solver itself reports completion and small algebraic residuals. Boundary value errors of zero do not imply PDE constraint satisfaction at boundary endpoints.

Next steps: independently recompute the fixed-band audit from the saved field arrays; inspect the one-sided derivative stencil and radial band sampling/weighting; compare truncation and modeling error without adjusting the convergence threshold post hoc. Keep the 65-grid candidate and artifacts frozen for review. Do not promote or tag the candidate as validated.
