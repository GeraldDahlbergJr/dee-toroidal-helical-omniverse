# Run 43 final numerical certification attempt

**Decision: NOT YET ESTABLISHED.** Independent audit agrees. Zero of the
60 scores has a justified total error within 10%, and zero of the 40 orders
has an unconditional certificate. This is missing proof, not demonstrated
failure of the underlying numerical solution.

Frozen parent: `c5580eea5b216052ba5e91c6725e0304a7852a99`.
Working branch: `run43-final-certification-20261107`. Execution date:
2026-11-08 UTC. All additions are in this directory; frozen evidence and
coefficients are unchanged. Three separate agents handled residuals,
nonlinear assembly, and independent review. No PDE solves were performed.

## Actual findings

* Every entry in both frozen SHA256 manifests matches. All three coefficient
  arrays have the expected finite binary64 shape.
* All 60 Gauss-6 and all 60 Gauss-8 archived scores reproduce bit for bit.
  A separate verifier replay also reproduces all 60 Gauss-8 scores exactly.
  These checks reproduce the same evaluator, not an independent enclosure.
* Exact rational sampled inverse upper bounds, rounded upward for display,
  are 1.323340785264946, 1.246889559229842, and 1.2364542730375343 for
  n=8,16,32. Their scope is the mathematical Gauss-6 Jacobian at the saved
  coefficients, under the frozen sampled-form derivation.
* The continuous denominator has new exact rational enclosures for all 60
  scores in the ideal protocol chart/spline space. Positivity and partition
  of unity bound psi; integrating J gives the geometric factor exactly.
  Machin's formula proves the pi bracket. The largest D_upper/D_lower is
  approximately 1.0903432728740536. Binary64 knot and evaluation discrepancies
  still need enclosure.
* All 40 conservative comparisons pass conditionally on valid 10% endpoint
  bounds. The minimum conditional order is 2.6643751474494595. Exact rational
  squared-ratio comparisons test p>=1.5; displayed logarithms are diagnostic.

## Pilot and obstruction

The n=32 first-outer-plane My score is 5.722436413449492e-7. Its exact
permitted total absolute error is the binary64 score divided by ten,
approximately **5.722436413449492e-8**. The q6/q8 gap is
3.6089171206474615e-14; it supplies no truncation certificate.

A rigorous, cancellation-destroying pointwise estimate yields the fixed-field
RMS range [0,115417926400], which is unusable. It excludes exact-assembly
solution uncertainty. The residual workstream derives a concrete sharper
test: bounding the sum of the two angular 16th-derivative suprema of the
weighted squared My integrand by 6.542024184414019e19 would bound Gauss-8
numerator truncation by 2.1551144541743275e-13. This is one illustrative
allocation only. The derivative bound has **not** been established.

The obstruction is preserving cancellation between spline second derivatives
and manufactured forcing while enclosing continuum integrals. PPoly conversion,
source and trigonometric evaluation, quadrature nodes/weights, derivative
arithmetic and reductions have no sharp outward reference enclosure here.
The next decisive residual test is a cellwise interval/Taylor enclosure of
that pilot numerator, including the quadrature remainder and all arithmetic
errors. A loose continuum range or additional empirical quadrature agreement
cannot substitute for it.

At n=32, no-solve Gauss-6 weak-defect dual-norm formulas evaluate in binary64
to 3.9976927784525283e-16 (scalar) and 1.4719210850224505e-14 (vector).
The direct q8-minus-q6 defect formulas are 7.236584956832922e-16 and
1.6120768955940396e-16. None is an outward-enclosed functional norm.
Exact-minus-Gauss assembly, absolute algebraic roundoff, Jacobian
perturbation and positive-ball nonlinear variation remain unbounded for
certification. No verified nonlinear solution neighborhood exists in this
evidence. The nonlinear report gives explicit contraction inequalities and
the additional surface second-derivative propagation requirements.

The decisive assembly test must enclose F and DF relative to Gauss-6 in
the stated dual/operator norms, establish a positive contraction ball, and
bound each strong residual's change on that ball. The total pilot allowance
must cover this propagated uncertainty as well as fixed-field integration.

## Evidence map

* `error_budget.json`: all 60 scores, exact allowances, denominator enclosures,
  separate unresolved arithmetic/assembly/solver/remainder entries, and all
  40 conditional comparisons. A null is unknown, never zero.
* `decision.json`, `decide.py`: fail-closed acceptance accounting; 0/60 and
  0/40 certified, 40/40 conditional. Machine bookkeeping cannot validate a
  future proof merely because someone fills in its fields.
* `residual/METHOD.md`: rigorous coarse pilot proof, sharp Gauss target, and
  q6/q8 reproduction artifacts and scripts.
* `nonlinear/REPORT.md`: rational inverse checks, empirical weak-defect
  diagnostics, and exact missing root/assembly/propagation obligations.
* `audit/REPORT.md`: independent mathematical and numerical review, rational
  coefficient-difference checks, denominator proof, and score replay.
* `REPRODUCE.md` and `SHA256SUMS`: bounded reproducibility and artifact seal.

The frozen verdict remains unchanged. Original-source production solves,
submission, baseline replacement, and fully validated solver status remain
outside this work.
