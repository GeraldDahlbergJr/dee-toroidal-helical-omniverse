# Run 43 cellwise interval/Taylor feasibility pilot (partial)

**Status: NOT ESTABLISHED.** The single n=32 first-outer-plane My total
error remains uncertified. No archived score or acceptance criterion changed.

## Implemented primitive

The bernstein.py module represents bivariate polynomials exactly as rational
monomial coefficients over a normalized angular cell. It simplifies algebraic
sums before range bounding, retaining exact polynomial cancellations.
Its conversion to Bernstein coefficients gives a rigorous convex-hull interval
for any polynomial on the cell. Binary64 inputs are embedded as exact rational
numbers, not interpreted as uncertain real-valued model parameters.

Unit tests include cancellation to zero, rational embedding, a conservative
Bernstein bound, and high-order polynomial differentiation.

## Exact boundary of this result

The Run 43 integrand is not a polynomial: the geometric frame involves
trigonometric and reciprocal functions; the manufactured field has
exponentials/trigonometric factors. The saved B-spline fields also require
exact piecewise coefficient extraction and normalized-cell conversion.
Those operations have not been incorporated into a rigorous Taylor model.
No 16th derivative supremum, denominator enclosure, quadrature error,
floating-point error, or propagated nonlinear solution error is obtained.

The first missing implementation obligation is a verified exact conversion
of the saved n=32 tensor-product spline restrictions and manufactured-source
terms to one common cellwise Taylor/interval expression with rigorous
trigonometric/exponential remainders. Naively bounding the components
separately would repeat the known cancellation-loss obstruction.

This implementation is a mathematical primitive, **not** evidence that
any of the 60 Run 43 scores satisfies the 10% error criterion. No broader
solver work or numerical simulations were attempted.
