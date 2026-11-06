# Run 43: certified element integration gate

Checkpoint: `ade8206d0f629fd8a73e6f66ef6490f0b1f5084d`.
One agent; one element; no mixed-operator assembly or physical solves.

**CERTIFIED scope:** enclosures of every entry of the full independent local
98-dimensional tensor mass matrix, and symmetry and positive definiteness of
that **exact** matrix, with a quantitative lower eigenvalue certificate.
The stored midpoint matrix is not certified positive definite. The enclosures
are deliberately coarse; this is a reproducible integration procedure, not an
assembled approximate Galerkin method or a coupled-solver certification.

## Exact integrands

Only the consolidated specification §§2,4,5 and conforming mixed proof §§1–4
are needed. No physical-source or production-solver definitions are needed
for this tensor mass. For an axis-aligned parameter brick with side lengths
`h_r,h_theta,h_phi`, set `d=h_r h_theta h_phi`, `q=1+r cos(theta)`.
The columns of `B=DF` are `h_r e_r, h_theta r e_theta, h_phi q e_phi`.
Consequently

    |det B| = d r q,
    g(f) = B^(-T) grad_xi(f)
         = e_r f_xi/h_r + e_theta f_eta/(r h_theta)
           + e_phi f_zeta/(q h_phi).

There is no Piola factor and no conformal `psi^6` volume weight. In particular,
the `1/r` and `1/q` factors inside the generators cannot be dropped. Let

    D_i(g) = e_i tensor g + g tensor e_i - (2/3)g_i I.

For reference monomials `f=xi^alpha`, `h=xi^beta`, with each exponent 0..2,
and enrichment functions `p,s` drawn from `1,F_x,F_y,F_z`, **all** tensor
mass entries are integrals over the reference cube of these three templates:

    d r q [2 delta_ij g(f)·g(h) + 2 g(f)_j g(h)_i
           - (4/3) g(f)_i g(h)_j],
    d r q [2 p(F) (E_beta g(f))_i],
    d r q [(E_beta:E_gamma) p(F) s(F)].

The STF Gram matrix is `diag(2,6,2,2,2)` in the specified E ordering.
The program uses all nine Cartesian entries through six symmetric components
with weights `(1,1,1,2,2,2)`, exactly equivalent to these formulas.
There are 81 D generators and 20 E generators. Only the three D generators
of the constant monomial are removed. Independence of the other 98 is
certified below, so this retains precisely the specified space on this cell.

For clarity, the derivative–derivative templates include geometric terms
proportional to `q/r`, `r/q`, `rq` and mixed inverse factors before contraction;
they are not polynomial degree-four integrands. There are no face integrals
in this tensor mass calculation. Its volume Jacobian does not replace the
specification's cofactor face Jacobian in any future boundary calculation.

## Selected method and exact-integration assumptions

The single selected integration method is validated interval subdivision.
The specification defines exact Lebesgue integrals, **not** mandatory symbolic
closed forms. We leave those exact targets and all equations unchanged and
return enclosures of them, rather than replacing them with a finite quadrature
bilinear form. No assumption is changed for that interpretation.

Analytic transformations remove some factors or integrate individual terms,
but this investigation does not establish closed-form antiderivatives for the
entire list, nor claim their impossibility. A symbolic-only implementation
requirement would need separate work: these enclosures are not symbolic exact
values. Ordinary Gaussian quadrature has not been used or claimed exact.

For a partition into `n^3` reference subcubes C, evaluate the **entire** reference
integrand `w T_i:T_j`, `w=d r q`, using Arb balls over each closed C. If `I_C`
contains its range, then

    H_ij in sum_C |C| I_C = [m_ij - R_ij, m_ij + R_ij],
    |H_ij - m_ij| <= R_ij.

This follows directly from monotonicity of integration; no sampling error
estimate is substituted for range inclusion. Arb outward bounds include
rational input rounding, transcendental rounding, products, sums and the
subcube-volume multiplication. No series is truncated, so there is no omitted
truncation term. The sole approximation is explicitly bounded by the stored
radius. Lower-triangle intervals are copied from the upper triangle.

`result.json` stores **every** upper-triangle midpoint and radius as exact
dyadics `[integer, exponent] = integer * 2^exponent`. Decimal display precision
therefore cannot silently reduce an error bound. Symmetry holds pointwise for
the exact integrand and is preserved identically by the stored enclosures.
Arbitrary independently chosen values from those intervals need not be
symmetric; they are not the claimed matrix.

All generators are analytic near this closed cell, and r and q stay bounded
away from zero. Natural interval extensions converge as the subcubes shrink
and working precision increases. This proves a convergent enclosure procedure
for the exact entries. The executable deliberately limits this bounded test
to 1..8 subintervals per axis and uses 128-bit integration precision. No
mesh-refinement family is generated.

## One admissible element

    A(xi,eta,zeta) = (1/10 + xi/100, eta/100, zeta/100),
    K = X(A((0,1)^3)).

The parameter side lengths are all 1/100. It lies in the original shell and
uses its exact curved geometry. It can be embedded in an admissible tensor
partition: use 20 radial intervals of length 1/100; for each angular circle,
use `[0,1/100]` and partition the remainder into 627 equal intervals. Their
length `(2*pi-1/100)/627` lies between 1/100 and 101/10000. Match both periodic
seams. Thus all cells have uniformly comparable side lengths and no hanging
faces. This describes admissibility only; no mesh is constructed or assembled.

## Positive definiteness: validated sample certificate and explicit eigenbound

This certificate concerns exact represented functions, not an inference from
positive diagonal entries or from the interval midpoint. Take 64 reference
points with each coordinate in `{1/5,2/5,3/5,4/5}`. Stack their six tensor
components in S, multiplying off-diagonal components by sqrt(2). Then

    A = S^T S,
    c^T A c = sum_s |T_c(x_s)|_F^2,    T_c = sum_i c_i T_i.

The program evaluates S at 1024 bits with outward bounds, forms an enclosing
Arb matrix A, and computes its validated inverse. A finite validated inverse
certifies nonsingularity of the exact sample Gram matrix. As that exact
matrix is positive semidefinite, it is positive definite. In particular, the
98 functions are independent. With the three discarded generators identically
zero, this also certifies the quotient of the entire supplied generator list.
No threshold-based numerical rank decision is used.

The sum of the absolute values of all inverse entries bounds its spectral
norm. Interval arithmetic certifies this sum is strictly less than `2^77`.
Therefore

    lambda_min(A) > 2^-77.

For any real coefficient vector c with Euclidean norm 1, at least one of the
64 sample points satisfies `|T_c(x_s)|_F > 2^-42 = alpha`, since
`alpha^2 <= 2^-77/64`.

Here are explicit conservative derivative and volume bounds transferring
this sample certificate to an integral certificate. On the reference cube,
`r >= 1/10`, `|q| >= 89/100`, `|partial f| <= 2`, and
`|partial partial f| <= 4` for every reference monomial. Frame components
have magnitude at most 1; their reference derivatives have magnitude at most
1/100. In the displayed g formula, coefficient magnitudes are bounded by
`100,1000,113`. Their reference derivatives are bounded by `0,100,2`.
The product rule gives each component derivative of g magnitude less than
8192, and each full Cartesian component derivative of D magnitude less than
32768. For the enrichment, `|F_i| <= 111/100`,
`|partial F_i| < 1/50` and `|E_beta,ij| <= 2` give much smaller bounds.
Cauchy–Schwarz over 98 generators, 9 components and 3 coordinate derivatives
therefore gives, conservatively,

    (sum_k |partial_k T_c|_F^2)^(1/2) <= L ||c||_2,    L=2^30.

Take `delta=2^-75`. A cube of half-side delta about any sample point lies
inside the reference cube, and `sqrt(3) L delta < alpha/2`. Throughout that
cube, `|T_c|_F >= alpha/2`. Finally,

    w = r q / 10^6 >= (1/10)(89/100)/10^6 > 2^-24.

Integration over that cube gives for every unit c

    c^T H_T c >= 2^-24 (2 delta)^3 alpha^2/4 = 2^-332.

Thus **lambda_min(H_T) >= 2^-332** in this specific 98-generator coefficient
basis. This is a deliberately conservative rigorous bound, not a condition
number estimate, a mesh-uniform eigenbound, or a bound for any global mixed
matrix. The inverse bound and the exponents are checked by the executable.

Entrywise bounds also control future numerical perturbations. With symmetric
midpoint M and radius matrix R, `|c^T(H_T-M)d| <= |c|^T R |d|`. In particular,
`||R||_infinity < 2^-332` would suffice to certify M positive definite and
supply a quantitative relative form error. **The preserved coarse enclosures
do not meet that condition.** We certify the exact mass target with the
sample certificate, not positive definiteness or norm equivalence of a
quadrature replacement. Global variational-error and face-compatibility proofs
would still be required before substituting a numerical form in the original
proof-only mixed method.

## Reproduction and tests

Tested environment: CPython 3.12.14; python-flint 0.8.0; bundled FLINT 3.3.1.
From the repository root:

```sh
python -m pip install -r research/run43_integration/requirements.txt
python research/run43_integration/certify.py --subdivisions 4 --output research/run43_integration/result.json
```

The JSON contains dependency versions, precision, all 4,851 independent entry
bounds, selected generators, the inverse/eigenvalue certificate, and passed
tests. In the preserved execution all entrywise absolute errors are strictly
less than `2^-1`; the individual exact dyadic radii are usually smaller and
are the authoritative entrywise bounds. Tests independently compare constant
STF diagonal entries with the
closed-form exact element volume

    V = (1/100)[(1/100)(r1^2-r0^2)/2
              + (r1^3-r0^3) sin(1/100)/3],
    r0=1/10, r1=11/100.

The E1, E2 and xy-pair mass entries enclose respectively `2V,6V,2V`; the
E1–E2 entry encloses zero. Further checks verify identical reflected
intervals, positive Jacobian and trace-free tensor values. The validated
sample inverse provides the nontrivial all-generator independence test.

The preserved run passed all assertions. A second execution into a temporary
file compared byte-for-byte equal with `cmp`. The result JSON SHA-256 is
`f9b1e47b35608f4cc2028e0909ebe4d4c36caf437e488233c3e52b9f53a7cd57`.

No production, original solver, source, audit or existing specification file
is modified. No scalar/displacement operators, trace constraints, coupled
loads, physical solves, nonlinear iterations, GMRES, refinement sweeps,
Actions or frozen-audit reruns are implemented or executed.

**One next action:** obtain authorization before mixed-operator assembly.
STOP.
