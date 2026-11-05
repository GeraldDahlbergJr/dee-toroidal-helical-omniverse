# Run 43 — far-field derivative design decision

Verified base: `4dd6277fb34f45dd9acd8622a724253a5c26b23b`.
One agent. Existing reports and archived certificates were read, not rerun.
This report uses symbolic reasoning only: no stencil, candidate space, numerical
experiment, physical solve, GMRES, grid sweep or GitHub Actions is constructed
or run.

**Centered far-field obstruction: CONDITIONAL.** The contradictions are
ESTABLISHED within their stated families; replacing every centered bulk
derivative is not proved mathematically necessary. **Recommended framework:
conforming mixed variational strain–displacement discretization with an exact
weak adjoint (variational summation by parts).** This selects a framework for
proof, not a successful replacement operator or an audit-passing method.

## What the existing contradictions actually require

| Evidence | Dependence and consequence |
| --- | --- |
| [Four-row derivation and certificate](run43_four_row_boundary/derivation.md) | Trapezoidal scalar SBP, affine exactness, collocated strain and fixed centered exterior force t=1/3 at the first two quadratic composition rows, then residuals 1/18 and -1/6. This rules out that family, not arbitrary derivatives. |
| [Six-row certificate](run43_six_row_controlled/certificate.json), [positive-norm derivation](run43_positive_boundary_norm/derivation.md), [variable-transition report and certificate](run43_variable_transition/report.md) | The weighted residual invariant is sum i W_ii R_i=-1/2. Allowing all skew couplings in the eight-node block, including two transition rows, does not remove it. Positive boundary weights also fail in their separately tested six-row family. The fixed nearest-neighbor centered outgoing interface, SBP and affine moments enter the proof. Changing the derivative/interface can invalidate this invariant; it does not guarantee feasibility. These are bounded-family certificates, not a classification of all SBP norms, supports or interiors. |
| [Original-adjoint report](RUN43_ADJOINT_CANDIDATE_20261105.md), [endpoint SBP derivation](run43_boundary_strain/derivation.md) | Original one-sided endpoint/centered columns and trapezoidal transpose give O(1/h) boundary errors; the minimal SBP endpoint replacement leaves O(1) error. These formulas depend on those particular closures, pairings and compositions. Changing the bulk alone need not repair them. |
| [Quadrature no-go](RUN43_QUADRATURE_NO_GO_20261105.md) | The original centered angular first derivative kills the checkerboard whereas the direct second derivative does not. Equal pointwise strain energy and unequal expanded work exclude every positive nodal reweighting of those fixed operators. A radial-only change does not resolve this angular mismatch. A redesigned strain/adjoint pair is outside that certificate. |
| [Product-rule obstruction](RUN43_PRODUCT_RULE_OBSTRUCTION_20261105.md) | Every linear derivation of the full finite nodal algebra with ordinary pointwise products is zero (apply Leibniz to nodal idempotents). This survives every nonzero replacement derivative in that algebra. The particular shifted centered identity does not survive unchanged. Exact ordinary Leibniz is neither available nor needed for a weak adjoint. |

The four-row t=1/3 kernel is an additional independent contradiction:
u=(0,1,-1,0,...) has Du=0, so angular-constant Cartesian W=u c has LW=0 on
the toroidal geometry despite zero radial trace. Its precise vector need not
remain a kernel of another derivative. What remains universally is

    ker(L0^T Ht L0) = ker(L0),    Ht>0, L0=Lh I.

Thus positive reweighting or a traction change cannot repair any retained zero
strain: EI=0 removes boundary work. Angular checkerboard derivatives alone do
not prove a full radial-Dirichlet strain kernel. Injectivity and uniform Korn
coercivity must be distinguished. Historical physical gate failures remain
failures; none is an impossibility theorem for another framework. See also
[feasibility report](RUN43_FEASIBILITY_ONLY_20261105.md) and
[completed closure report](RUN43_COMPATIBLE_CLOSURE_20261103.md).

## One justified framework, and why

Use conforming displacement spaces Vh contained in H1 on the actual shell,
periodic in the angles, with Vh,0 having zero trace at both radial faces. Use an
STF tensor space Th containing the actual Cartesian strain L(Vh), with exact
volume pairings. The first derivative is the weak derivative of the conforming
field, not a collocated centered matrix. Introduce strain sigma through

    (sigma,tau) = (L W,tau)       for every tau in Th,
    -(sigma,L V)/2 = (f,V)      for every V in Vh,0,
    f = 8 pi G psi^6 S.

Containment makes sigma=LW exactly; this is an auxiliary variable, not extra
physical data. Eliminating it gives precisely the original weak momentum
equation restricted to Vh,0. This is a framework specification, not a chosen
element, degree, mesh, quadrature or candidate implementation.

There is a concrete kernel argument beyond merely hoping staggering works.
For zero-trace Cartesian fields in three Euclidean dimensions, integration by
parts gives

    ||L V||_L2^2 = 2||grad V||_L2^2 + (2/3)||div V||_L2^2.

Boundary terms vanish; periodic seams contribute none. Density extends the
identity to H1 zero-trace fields. Poincare then bounds the H1 norm by the strain
norm with a domain constant independent of h. Every conforming subspace with
exact strain and integration inherits this bound. Therefore a nonzero conforming
representation of the old mode cannot have zero strain. Projecting strain into
an insufficient tensor space or underintegrating it would lose this argument.
This proves an abstract route around the known kernel, not stability of a
particular assembled toroidal operator or of the nonlinear coupled system.

A staggered pair can also satisfy rectangular SBP identities, but placement
alone proves neither kernel removal nor consistency. Existing evidence selects
no staggered stencil. The conforming mixed framework is recommended because
its containment and integration hypotheses supply the explicit kernel proof
above while avoiding the failed same-matrix D-squared restrictions.

## Minimum obligations before implementation

1. **Polynomial and approximation conditions.** State actual spaces, mappings
   and pairings. For derivative maps between representations Pv and Pg require
   Dh Pv p=Pg p' for reproduced polynomials, at least p=1,x; for Cartesian
   gradients require Gj Pv p=Pg partial_j p. State the reproduced degree and
   prove approximation estimates. For a strong nodal normal composition
   C=Divh Lh reduced to its scalar principal factor, require C 1=0, C x=0,
   C x^2=2 at every unconstrained boundary and transition row. For order p
   with fixed scaled coefficients require C x^k=k(k-1)x^(k-2) through k=p+1;
   second order includes the cubic moment. These tests use local smooth
   Dirichlet extensions, not inadmissible global polynomial boundary values.
   A weak formulation instead needs the exact polynomial weak residual
   -(L p,L v)/2=(div L p,v)-<v,(L p)n> and consistent loads. Weak exactness
   alone does not imply these strong boundary moments.
2. **Adjoint and geometry.** Prove positive volume and surface pairings and
   the full Cartesian STF Green identity
   Hv Divh + Lh^T Ht/2 = Ev^T Hb N on STF tensors (all nine entries counted).
   For a rectangular scalar pair this reads
   Hv Ddiv + Dgrad^T Ht = Ev^T Hb Et, with oriented flux trace Et.
   Then M A_D=-(Lh I)^T Ht(Lh I)/2. Use J=r(1+r cos theta) and the actual
   geometry and normals; scalar SBP alone is insufficient. Exact integration
   establishes the variational version; any computable quadrature must prove
   the same identity and uniform norm equivalence or a controlled consistency
   error, without assuming a nodal product rule.
3. **Kernel removal.** Prove strain containment, faithful displacement
   representation, ker(Lh I)={0}, and a mesh-independent Korn bound in the
   stated H1-equivalent norm. If strain is projected, prove the corresponding
   uniform lower bound for that projection; detection of only the old mode
   is insufficient. Prove stability of the auxiliary elimination as well.
4. **Boundary consistency and unchanged acceptance.** Prove approximation to
   the original Dirichlet problem and consistent boundary reactions/traction.
   Before claiming a strong boundary replacement, prove vanishing truncation
   error at all affected layers, including both faces and transitions, with
   variable geometry. Weak convergence alone cannot meet that claim. Specify
   consistent reconstruction for the unchanged frozen audit and Hamiltonian
   strain coupling; prove relevant error control without redefining either.
   No claim that the twelve-check gate will pass follows from this report.

## Preservation and stop

The original Einstein–DEE equations, action, sources and Hamiltonian remain
those in the repository: Delta psi+psi^-7|LW|^2/8+2 pi G psi^5 E(psi)=0 and
div LW=8 pi G psi^6 S, with the original Cartesian STF L. Retain G=.001,
a=1/20, b=1/4, psi=1 and W=0 radially, and angular periodicity. No zero
traction, derivative condition or source compatibility constraint is added.
Production, frozen audit and historical evidence are unchanged. Only this
report is added from the verified base.

One next necessary step: separately authorize a proof-only specification of
conforming displacement/tensor spaces and integration satisfying the four
obligations above. No new candidate or physical solve is authorized here. STOP.
