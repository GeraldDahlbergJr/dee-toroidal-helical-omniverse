# Run 43 targeted evaluator qualification repair — FAIL

Parent: `8d9070429c5f9f77c70f89858d3f6dcb3ec4fc1a`; branch `validation/run43-final-bounded`. One agent, one targeted repair, zero DEE solves. Numerical work stopped at qualification failure. No baseline residual reevaluation or additional repair cycle.

| Evidence | Status | Basis |
|---|---|---|
| Radial coefficient construction repair | PASS | Same nine-point Lagrange derivative; exact rational uniform offsets and weights, converted once to float |
| Constant annihilation, n=8,16,32,64 | PASS | Both first and second derivatives: 8/8 under unchanged gamma dot-product bound |
| Polynomial derivatives | PASS | Degrees 0–8, both derivatives, all four grids: 72/72 under coefficient-conversion plus application roundoff bound |
| Manufactured H, Mx, My, Mz qualification | FAIL | Every component fails p>=1.5 on 16->32 and 32->64 |
| Accepted inputs and equations unchanged | PASS | All 31 preexisting manifest hashes match; no solver or archive edits |
| Baseline continuum residuals, uncertainty, boundary checks and conservative rates | NOT YET ESTABLISHED | Conditional evaluation withheld because qualification failed; historical diagnostics remain unqualified |
| PRD readiness | NOT YET ESTABLISHED | Do not claim independent continuum validation or recommend PRD submission on this evidence |

## Diagnosis and remaining qualification defect

The old construction formed floating-point monomial coefficients from nearly integral stencil offsets and normalized each basis by evaluating that polynomial at its node. Cancellation in coefficient construction and normalization broke the exact Lagrange moment identities. Exact rational weights for the same uniform stencil remove that construction error; the stored matrices now pass the original constant test without row-sum adjustment or tolerance relaxation.

The remaining precise qualification defect is that the retained Cartesian cubic fixture is analytically resolved by degree-eight radial interpolation and the angular Fourier grid: its scalar pullback has radial degree at most three and angular frequencies at most three; its Cartesian strain has degree two. Consequently it provides no nonzero truncation-error sequence on these grids from which p>=1.5 can be demonstrated. Floating-point differentiation of the sampled fields still produces nonzero errors that grow with refinement, including direct second derivatives of psi near one. Polynomial checks establish the matrix moments, but do not bound accumulated sampling, metric, Fourier and composed-strain roundoff in the complete residual. No full residual roundoff bound was previously established, so no roundoff exemption or manufactured absolute tolerance is introduced. These results do not prove another physical-equation defect, nor a failure of the accepted baseline. A qualified continuum assessment remains unavailable.

## Manufactured errors

| n | max H error | max Mx error | max My error | max Mz error |
|---|---|---|---|---|
| 8 | 3.070538268e-10 | 8.428240744e-14 | 8.127266221e-14 | 1.099120794e-14 |
| 16 | 1.863187826e-09 | 3.659312436e-13 | 4.065237730e-13 | 9.593888184e-14 |
| 32 | 8.575841148e-09 | 2.291103071e-12 | 2.830254260e-12 | 8.404531411e-13 |
| 64 | 4.269892973e-08 | 1.004234022e-11 | 1.017137546e-11 | 2.179055764e-12 |

| pair | p H | p Mx | p My | p Mz |
|---|---|---|---|---|
| [16, 32] | -2.202505 | -2.646398 | -2.799520 | -3.130980 |
| [32, 64] | -2.315850 | -2.131981 | -1.845511 | -1.374464 |

## Reproduction and scope

Run `OPENBLAS_NUM_THREADS=1 python research/run43_final_validation/test_manufactured.py` from the repository root. Expected exit: 1. It writes only `repair_qualification.json`; redirect stdout to `repair_qualification.log` to retain the process evidence. It computes the manufactured residual directly and never loads saved fields.

Polynomial test inputs s=i/n and s**k are exact dyadic numbers on these grids; their physical derivatives are exactly representable here. The bound is (eps/(1-eps)+gamma_(n+1))*abs(D)@abs(v), covering conversion of exact rational coefficients and stored-matrix multiplication. This is an arithmetic bound, not a physical tolerance.

Historical REPORT.md, results.json, qualification.json and their original checksum manifest are preserved as prior failed evidence; they are not outputs of this repair. Do not run the historical renderer to describe this attempt. New evidence is in this report, repair_qualification.json/log and repair_input_integrity.json. The evaluator main routine remains the historical unqualified diagnostic path and was not executed.

PRD recommendation: NOT READY on the available independent continuum evidence. STOP; await further authorization.
