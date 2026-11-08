# Run 43 final manufactured uncertainty audit

Evidence parent: b75afe94079f61d890712e06208d5226f312897e.

**Verdict: NOT YET ESTABLISHED.** All 60 score budgets are enumerated.
The archived evidence supplies no complete total-error enclosure for any
score. This audit performed no PDE solves and did not modify coefficients.

| Requirement | Decision |
|---|---|
| Exact Gauss-6 coupled Jacobian inverse at the archived fields | PASS. Uniform bound below 1.33 in the max of scalar/vector gradient norms, with the corresponding max dual input norm. |
| Fixed-field quadrature truncation | NOT YET ESTABLISHED. The measured q6/q8 proxy is empirical. |
| Floating-point evaluation and derivative arithmetic | NOT YET ESTABLISHED. No outward-rounded strong evaluator reference enclosure. |
| Propagated assembly and solver error | NOT YET ESTABLISHED. No exact-assembly weak defect or verified root neighborhood. |
| 10% total error | NOT YET ESTABLISHED for all 60 scores. |
| Conservative convergence | Conditional PASS for all 40 comparisons. With valid 10% enclosures, minimum orders would be 2.6740485344 and 2.6643751474. Unconditional orders remain unestablished. |

## What was verified

The inverse bound follows from the existing certified scalar/vector
coercivities a,v, cross bounds c,d, and Schur bound s. Block elimination
gives B_scalar=(1+c/v)/(a(1-s)) and B_vector=(1+d B_scalar)/v.
The largest is below 1.33 on each grid, using exact rational arithmetic on
the archived outward-rounded inputs. This applies at the saved coefficients
to the mathematical Gauss-6 Jacobian.

The exact permitted absolute error is 0.1 times each archived binary64
score. For any comparison whose endpoint scores both meet that allowance,
the guaranteed refinement ratio is 0.9 E_coarse/(1.1 E_fine). The script
tests its square against 8, avoiding rounding in the p>=1.5 comparison.
Displayed logarithms are approximate. An independent review verified
both minima directly from the archived scores. A separate review derived
the inverse and defect propagation requirements.

## Missing evidence

For the stored fields, enclose the squared residual integral N and the
volume denominator D, including strong derivative arithmetic and quadrature
truncation. For example, an independently outward-rounded reference integral
plus a cellwise analytic Gauss remainder bound would permit direct comparison
to the archived score. Relative N and D bounds of 15% and 3%, if established,
give RMS factors [0.9084287,1.0888375] before adding any propagated error.
At n=32 the tightest absolute total allowance is 5.7224364e-8 (first outer
plane My). No such enclosures were computed in this audit.

For the separate exact-assembly solution target, bound the weak defect in
the gradient dual norm, including the absolute Gauss-6 algebraic defect and
exact-minus-Gauss-6 assembly difference. Enclose the assembly Jacobian
perturbation and nonlinear remainder in a positive neighborhood. Then bound
the resulting change in each strong residual norm, including surface
estimates on the first interior planes. The inverse bound at the saved
coefficients does not by itself supply those quantities. The relative
algebraic residual and strong-rule sensitivity cannot replace them.

## Reproduce

From the repository root:

    python research/run43_final_gate/audit.py research/run43_bspline_completion > /tmp/run43-final-gate.json
    cmp /tmp/run43-final-gate.json research/run43_final_gate/audit.json
    (cd research/run43_final_gate && sha256sum -c SHA256SUMS)

Original-source execution and manuscript preparation remain gated.
