# Nonlinear and assembly uncertainty: NOT YET ESTABLISHED

Checkpoint: `c5580eea5b216052ba5e91c6725e0304a7852a99`.
No PDE solves or coefficient changes were performed. The sampled conditioning
calculation reproduces the archived JSON exactly. `certify.py` independently
forms the Schur product from the outward-rounded block coercivity/cross bounds
using exact rational arithmetic, rather than assuming the archived Schur entry.
The resulting upward-rounded max-product-gradient inverse bounds are:

| n | sampled inverse upper |
|---|---:|
| 8 | 1.323340785264946 |
| 16 | 1.246889559229842 |
| 32 | 1.2364542730375343 |

These bounds apply to the mathematical Gauss-6 Jacobian at the archived fields,
conditional on the analytic sampled-form derivation in the frozen CONDITIONING.md.
They do not bound the inverse of the exactly integrated Jacobian.
`bounds.json` records the exact rational bounds and explicitly null missing quantities.

## Absolute defect diagnostic, with scope limitation

`weak_defect.py` evaluates the weak scalar/vector defects at archived fields.
It never calls `solve`, `linear`, CG, or another PDE solver. Its binary64 results
are diagnostics; neither a root radius nor a rigorous solver-error bound follows.
The output explicitly states `rigorous: false`.

Let H denote the unweighted parameter tensor stiffness with coefficients
(1, 1/.15^2, 1). On the chart, J >= .0375,
J/r^2=q/r >=3, and J/q^2=r/q >=.04. Hence the physical scalar gradient Gram
matrix K satisfies K >= .0375 H, both for exact integration and mathematical
Gauss-6 (H is exactly integrated). The archived `pre` represents H^-1/.15,
so sqrt(4 r^T pre(r)) is a valid *formula* for an upper bound on a scalar
dual norm in exact arithmetic; vector component bounds combine in quadrature.
The experiment's floating-point eigenvectors, generalized eigenvalues, basis
and source evaluations, sums, and residual cancellations have not been
outward-enclosed. Accordingly the computed formula is not itself a certified
upper bound. It is useful for deciding where to invest interval arithmetic.

At n=32 Gauss-6, the diagnostic formula gives scalar 3.9976927784525283e-16
and vector 1.4719210850224505e-14. The Euclidean coefficient defect norms
are respectively 2.054060116666164e-17 and 4.9763212342884445e-15.
At Gauss-8 the corresponding formulas are 7.311859795267213e-16 and
1.471954762490448e-14. Direct Gauss-8-minus-Gauss-6 weak coefficient
defects give diagnostic dual formulas 7.236584956832922e-16 (scalar)
and 1.6120768955940396e-16 (vector). These are rule-difference diagnostics,
not upper bounds for exact-minus-Gauss-6 assembly.
The cancellation scale makes plain why small reported residuals cannot stand
in for an absolute assembly roundoff bound.

## Exact requirements for a verified neighborhood

Fix the max-product gradient norm X and its corresponding max dual norm.
Write F for the exact weak map, Q for mathematical Gauss-6 assembly,
x0 for the archived coefficients, and B for the sampled inverse upper above.
A usable proof would enclose

* eta >= ||F(x0)||_(X*) <= ||Q(x0)||_(X*) + ||F(x0)-Q(x0)||_(X*),
* delta >= ||DF(x0)-DQ(x0)||_(X -> X*),
* L such that ||DF(x)-DF(x0)|| <= L ||x-x0|| on a positive radius rho ball.

For A0=DQ(x0), the map T(x)=x-A0^-1 F(x) is a contraction on that ball if
k=B(delta+L rho)<1 and B eta+k rho<=rho. These inequalities give a verified
solution and local uniqueness. Every inequality needs rigorous upper bounds,
and positivity must hold throughout the ball. A sampled inverse at x0 alone
supplies none of eta, delta, L, or rho.

For each strong residual score, including the two first-interior surfaces,
a further trace/inverse estimate K_j is required with
||R_j(x)-R_j(x0)||_(L2 normalized) <= K_j ||x-x0||_X.
The resulting propagated uncertainty K_j rho must enter that score's budget.
A volume energy estimate alone does not control surface second derivatives
without these finite-dimensional constants.

## Precise obstruction and next decisive test

No outward-enclosed absolute Gauss-6 weak defect, exact-minus-Gauss assembly
functional norm, assembly Jacobian perturbation, positive-ball nonlinear
Lipschitz bound, or strong surface propagation constants are available.
The diagnostic binary64 residuals cannot fill these fields. Therefore no
verified nonlinear solution neighborhood is established and solution
uncertainty remains unbounded in the certification budget.

The next decisive assembly test is an outward-rounded cellwise weak-integrand
reference enclosure (including the manufactured source) with a rigorous
quadrature remainder and a certified conversion of coefficient defects to
X* norms. It must include DF as well as F. This test should first target
n=32 and then be coupled to surface inverse constants for the first-outer My
score; its entire final allowance is only approximately 5.7224364e-8.
Simply increasing Gauss order does not close any of these proof obligations.

## Reproduce

From repository root:

```sh
python research/run43_final_certification/nonlinear/certify.py > /tmp/run43-nonlinear.json
cmp /tmp/run43-nonlinear.json research/run43_final_certification/nonlinear/bounds.json
OPENBLAS_NUM_THREADS=1 timeout 240 python research/run43_final_certification/nonlinear/weak_defect.py 32
OPENBLAS_NUM_THREADS=1 timeout 240 python research/run43_final_certification/nonlinear/weak_defect.py 32 8
OPENBLAS_NUM_THREADS=1 timeout 240 python research/run43_final_certification/nonlinear/weak_defect.py 32 6 --compare
```

The optional last argument chooses diagnostic quadrature (default 6).
Timing is not a reproducibility requirement. NumPy 2.3.5, SciPy 1.17.0,
and SymPy 1.14.0 match the archived environment. No empirical diagnostic is
used as an acceptance bound. Frozen criteria remain unchanged.
