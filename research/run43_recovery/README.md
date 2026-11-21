# Run 43 mixed-flux + SBP–SAT recovery experiment

**Priority:** first numerical recovery experiment. **Separate next physics extension:** Wendelstein 7-X microwave–plasma pilot (`research/run43-w7x-flux-potential-pilot-20261120`).

## Guardrails
- Original frozen Run 43 **12/12** certified baseline remains unchanged.
- Experimental higher-order variation was previously reported as **9/12**. This experiment does **not** claim to fix those three checks.
- Branch originates from `research/run43-interval-taylor-pilot-20261111`; no changes to the original DEE solver.
- This is a minimal **synthetic scalar** verification problem, not full coupled Einstein–DEE evolution.

## Mathematical setup
Solve `-u'' + u = f` on `[0,1]` with manufactured solution
`u(x) = 1 + x + sin(pi*x)`, `f(x) = pi^2*sin(pi*x) + u(x)`,
and endpoint Dirichlet values `g=(1,2)`. Auxiliary flux `j=D u`.

The diagonal-norm second-order-interior SBP first-derivative satisfies
`H D + D^T H = B`. The eliminated mixed-flux weak/SAT operator is
`L = D^T H D - B D - D^T B + H + W`, where `W` is the boundary
penalty matrix (`12/h` at both endpoints). Its load is
`b = H f - D^T B g + W g`.

Here the two symmetric boundary-consistency terms are **essential**;
a penalty-only variant demonstrated approximately first-order behavior
and is **not** used in the checked-in test.

## Reproduce
`python research/run43_recovery/mixed_flux_sbp_sat.py`

With NumPy and BLAS available, check n=16,32,64,128, the discrete SBP
identity, positivity of the weak operator, solution residual and
manufactured solution convergence `p >= 1.5`.

## Initial local results (2026-11-20)
| n | L2(H) solution error |
|---|---|
| 16 | 0.0083085043 |
| 32 | 0.0020662967 |
| 64 | 0.0005159102 |
| 128 | 0.0001289364 |

Observed p = 2.00754, 2.00186, 2.00046; SBP identity zero to
floating-point precision. These are **pilot** results only.

## Next real Run 43 gate
Reproduce exact 9/12 experimental failures; identify failing tests and
boundary operators. Map the genuine DEE radial flux and physical boundary
conditions into a compatible SBP–SAT discretization and prove discrete
stability/coercivity before re-running all twelve checks at the frozen
8/16/32/64 resolution hierarchy. Do not loosen acceptance tolerances.

After this numerical repair, the W7-X microwave–plasma extension must
introduce a physically justified electromagnetic response model; the
synthetic magnetic potential pilot alone is not a plasma simulation.
