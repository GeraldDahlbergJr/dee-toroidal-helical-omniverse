# Independent audit: NOT YET ESTABLISHED

Frozen source: `c5580eea5b216052ba5e91c6725e0304a7852a99`.
No PDE solves were performed. This audit independently reads archived binary64
coefficients and scores; it does not rely on the other workstreams' verdicts.

`verify.py` checks both archived SHA256 manifests with actual file reads. Every
listed file matches, including the two historical log files. All three archives
have the declared finite binary64 shape. Exact rational coefficient subtraction
independently verifies the archived outward coefficient, derivative, coercivity,
cross-coupling and Schur inputs. Exact rational block elimination reproduces the
sampled inverse caps below 1.33. These are bounds for the mathematical Gauss-6
Jacobian at the saved fields, conditional on the frozen analytic sampled-form
proof. They do not enclose exact assembly or a nonlinear solution.

All 40 conservative convergence comparisons pass **conditionally** on valid
10% endpoint error enclosures. Acceptance uses the exact rational test
`(9 E_coarse/(11 E_fine))^2 >= 8`, not rounded logarithms. The minimum displayed
conditional order is 2.6643751474494595. None has an unconditional certified
order because the endpoint enclosures are missing. The tightest score allowance
is first-outer-plane My at n=32, 5.722436413449492e-8.

`denominator.py` independently verifies all 60 rational denominator enclosures.
It proves the supplied pi bracket using Machin's identity with exact rational
50-term alternating arctangent series. Positivity and spline partition of unity
give bounds on psi; integrating positive J against those bounds gives the
reported enclosures. The cos(theta) contribution integrates to zero. These
bounds address the protocol's ideal exact decimal chart and ideal spline
space. On first-interior surfaces the protocol uses J angular weights, as the
budget does; it does not request geometric surface area. No bound on discrepancy
between the ideal basis and binary64 PPoly evaluation follows from this check.

Independent q8 replay reproduces all 60 scores bit for bit (maximum relative
difference zero). It invokes the archived evaluator separately from Agent 1,
so it verifies execution and archive reproducibility, not an independent
mathematical implementation or an enclosure.

The residual pilot's Gauss-8 constant `(8!)^4/(17 (16!)^3)` and tensor composite
remainder form `L^2 C h^16 (M_theta+M_phi)` are valid for cellwise C16 numerator
integrands, with M bounds for the pure 16th derivatives. The displayed ceiling
is a sufficient target for its proposed numerator allocation, not an established
derivative bound. Its lower denominator estimate is conservative and the n=32
scalar coefficients are positive. Neither a numerator enclosure nor a complete
RMS error enclosure has been produced.

The nonlinear report correctly distinguishes empirical weak-defect diagnostics
from rigorous dual-norm bounds. Its contraction conditions require bounded exact
weak defect eta, Jacobian assembly perturbation delta, and nonlinear variation
L on a positive ball. Those bounds, positivity on the ball, and propagation of
solution uncertainty to second-derivative surface residuals are absent. A small
binary64 residual or sampled inverse cap does not supply them.

The additional coarse continuum pilot was independently rerun and reviewed.
The spline derivative recurrence gives first derivatives <=8/h0 and second
derivatives <=64/h0^2; the conservative frame and source constants imply the
fixed-field RMS enclosure [0,115417926400]. This is valid but unusably wide,
and is not a rigorous bound on the exact-assembly solution score.

The acceptance checker was reviewed: on the present null budgets it returns
0/60 certified scores, 0/40 certified comparisons, and all 40 conditional tests
passing. Its mechanically checked nonnull fields cannot replace independent
proof review in a future run.

The 60-row error budget correctly uses null to mean unestablished, not zero.
The exact-domain qualifier and missing arithmetic/assembly terms preclude any
10% claim. **Independent decision: NOT YET ESTABLISHED**, with unchanged criteria.
The next decisive test is a rigorously outward-enclosed n=32 first-outer My
squared-numerator integral including cellwise quadrature remainder and basis,
source, geometry and arithmetic errors; the separate exact-assembly root and
surface-propagation obligations must also be closed before certification.

Reproduce from repository root:

```sh
OPENBLAS_NUM_THREADS=1 python research/run43_final_certification/audit/verify.py > /tmp/independent.json
cmp /tmp/independent.json research/run43_final_certification/audit/independent.json
python research/run43_final_certification/audit/denominator.py > /tmp/denominator.json
cmp /tmp/denominator.json research/run43_final_certification/audit/denominator.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 240 python research/run43_final_certification/audit/replay.py > /tmp/replay.json
cmp /tmp/replay.json research/run43_final_certification/audit/replay.json
```

The replay invokes the archived strong evaluator independently of the residual
workstream. Its agreement is a reproducibility check, never an error certificate.
