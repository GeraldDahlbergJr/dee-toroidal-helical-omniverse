# Exact sampled Galerkin Jacobian certificate

This certificate addresses the mathematical six-point Gauss-per-span
Jacobian at the archived spline coefficients. It does not certify a
floating-point matvec error, a nonlinear root enclosure, an assembly
quadrature error, or an evaluator uncertainty bound.

## Principal sampled forms

Write Q for tensor Gauss-6 on the parameter mesh, J=rq for the volume
Jacobian, and use the energy norm ||v||_g²=Q(J|grad v|²).
The scalar stiffness is its squared norm. The Cartesian vector weak form
is a_Q(W,W)=Q(J|LW|²/2). Pointwise,

    |LW|²/2 = |grad W|² + (div W)²/3 + N(W),
    N(W) = tr((grad W)^2) − (div W)^2.

The exact physical integral of N is zero by Cartesian integration by
parts, using zero radial traces and the periodic seams. The sampled
integral need not vanish. The following bound is for the actual Gauss
rule rather than the exact variational form.

Put g_a=partial_a W and l=(1,r,q). The coefficient of the antisymmetric
products in J N is a cofactor of the reciprocal chart frame. Explicitly
J N is a signed sum of three terms
2 partial_c X dot (g_a cross g_b), one per unordered pair a,b, with c its
complement. The covariant frame vectors are e_r, r e_theta and q e_phi.
All denominators have canceled.

On each angular span approximate sin and cos by their midpoint Taylor
cubics. For n>=8, each scalar remainder is bounded by
(pi/8)^4/24<.001. The product of two approximants in different angles has
absolute error at most .002001. For e_r and e_theta the resulting
Euclidean error is less than .0031. For q e_phi, each nonzero component
error is bounded by .001 max(q)+.25*.001*1.001, and the Euclidean error
is also less than .0031 max(q), since q>=.75. Thus every covariant vector
approximation has error at most .0031 max_cell(l_c). Its polynomial has
degree at most 3 in each angle and at most 1 in radius. The products of
spline derivatives are of degree at most 8 in each variable, so their
product with the approximant is integrated exactly by Gauss-6.

Take w_a=min_cell(J/l_a²). Since
max_cell(l_c)/sqrt(w_a w_b) <= rmax qmax/(rmin qmin) <= 2.5
on every cell, the pointwise error in the null term is bounded by

    .0031*2.5 sum_{a<b}(w_a |g_a|²+w_b |g_b|²)
      = .0155 sum_a w_a |g_a|².

The polynomial at right is integrated exactly by Gauss-6, and its sampled
norm is bounded by .0155 Q(J|grad W|²). This bounds each of the exact and
sampled remainder integrals. Consequently
a_Q(W,W) >= .969 ||W||_g², with no unexpected sampled vector kernel.

For either scalar or vector homogeneous radial traces, unweighted radial
mass and stiffness are exactly integrated. The one-dimensional Dirichlet
Poincare inequality and the chart bounds imply
Q(J|v|²) <= (Jmax/Jmin)(.2/pi)^2 ||v||_g² < .034 ||v||_g².
This also establishes that the scalar gradient seminorm is a norm.

## Coupled Schur estimate

At the stored coefficients, let A be the scalar Newton diagonal block,
V the vector stiffness, and C,D the scalar/vector and vector/scalar
off-diagonal blocks. Since the strain contribution to the scalar Newton
potential is nonpositive, an upper bound on that potential is

    P = 2 pi G (5*.0127638*psi_max^4 + .068614375*(400+9/.75^2)).

It follows that A has lower bound alpha=1-.034 P in gradient norm.
For the off-diagonal blocks, use |L delta W|<=2|grad delta W| and
||S||_infinity < .62. The stated Poincare bound and Cauchy–Schwarz give

    ||C||_(g→g*) <= .185 sup|LW| / (2 psi_min^7),
    ||D||_(g→g*) <= .034*48*pi*G*psi_max^5*.62.

Therefore the induced gradient-norm bound on the scalar Schur correction
A^-1 C V^-1 D is ||C|| ||D||/(alpha*.969). A bound strictly below one
certifies invertibility by the Neumann series.

The script encloses psi by the nonnegative partition of unity (including
the omitted zero endpoint coefficients), and bounds each parametric
derivative component by 4/h times the maximum adjacent coefficient
difference. This follows from the standard derivative coefficient
identity and the degree-three nonnegative partition of unity. Angular
coefficient differences wrap periodically. The Cartesian gradient
Frobenius norm equals its orthonormal chart-frame value; the script
uses the sum of entrywise bounds. Rational arithmetic and conservative
pi bounds are used after converting stored IEEE coefficients and
outward-expanded floating-point coefficient differences.

Run from the evidence directory:

    python conditioning.py coefficients_8.npz coefficients_16.npz coefficients_32.npz

This tests only archives available after the authorized run. It performs
no solves. Rounded output bounds are outward-expanded.
