# Run 43 — one conforming mixed strain–displacement specification

Starting verified commit: `5c72071c0558fd54e168e3f062de6bad27aca567`.
Research branch: `research/run43-conforming-mixed-proof`.
One agent; one design; analytical derivation and existing evidence only.

**Decision:** the spaces below have exact strain containment, an exact
variational adjoint, and a mesh-independent Dirichlet Korn bound. They define
a consistent conforming discretization of the original momentum equation.
Compatibility with the unchanged nonlinear Hamiltonian and frozen twelve-check
strong-residual audit is **UNRESOLVED**. The precise missing implication is
identified in Section 7. This is the stopping point, not an audit-passing method
or permission to implement or solve.

The starting point is the [far-field design decision](RUN43_FAR_FIELD_DESIGN_DECISION_20261105.md),
especially its four minimum obligations. Its bounded-family contradictions and
historical failures are retained as evidence; none is rerun or reinterpreted as
a universal impossibility theorem. No production, source, audit, or historical
file is changed by this work.

## 1. Original domain and differential operator

Use Cartesian vector components throughout. With a=1/20, b=1/4 and R=1,

    X(r,theta,phi) = ((1+r cos theta) cos phi,
                      (1+r cos theta) sin phi, r sin theta),
    (r,theta,phi) in (a,b) x (real numbers / 2pi Z)^2,
    q = 1+r cos theta,  J = r q > 0.

The physical domain Omega is this exact toroidal shell, without faceting or
polynomial approximation of its boundary. Write er, etheta, ephi for the
orthonormal coordinate directions. Then

    grad u = er u_r + etheta u_theta/r + ephi u_phi/q,
    dV = r q dr dtheta dphi,
    n_a = -er, n_b = +er,  dS|_(r=c) = c q(c,theta) dtheta dphi.

All geometric coefficients and their derivatives are the actual analytic ones.
The map is nonsingular: r>=1/20 and q>=3/4. Angular seams are identified; they
are not physical boundaries. Define

    (L W)_ij = partial_j W_i + partial_i W_j
               - (2/3) delta_ij partial_k W_k,
    (div T)_i = partial_j T_ij,
    f(psi) = 8 pi G psi^6 S,   G=0.001.

The original equations are

    div L W = f(psi),
    Delta psi + psi^-7 |L W|^2/8 + 2 pi G psi^5 E(psi) = 0,
    W=0, psi=1 at r=a,b; angular periodicity.

The source S and physical E(psi) are exactly those in
[the production formulation](../reproduce/solve_toroidal_helical_coupled_constraints.py)
and [the original source](../reproduce/validate_toroidal_helical_dee_source.py).
There is no traction condition, extra source restriction, or change in coupling.

## 2. One concrete mapped element pair — SPECIFIED

Use conforming tensor-product partitions of the parameter shell, periodic in
theta and phi, with matching opposite seam faces. Each parameter brick is the
affine image A_K of the reference cube (0,1)^3. Require uniformly comparable
parameter side lengths and bounded aspect ratios, and refine with maximal side
length h tending to zero. This is sufficient for uniformly regular physical
curved elements K=F_K((0,1)^3), F_K=X composed with A_K. Uniform regularity follows
from the fixed smooth nonsingular geometry and the parameter mesh bounds.
No hanging faces are included in this specification.

Let Q2 denote reference polynomials of degree at most two in **each** variable,
and Pm denote physical Cartesian polynomials of total degree at most m.
The scalar local displacement space is the following explicitly enriched space:

    U(K) = { qhat composed with F_K^-1 : qhat in Q2 }
           + { p restricted to K : p in P2(R^3) }.

Thus its reference pullback is Q2 + span{1, F_i, F_i F_j : 1<=i<=j<=3}.
Redundant generators are identified as functions, not treated as independent
unknowns. The polynomial degrees are reference tensor degree 2 and physical
Cartesian enrichment degree 2. The trigonometric geometry is exact, of no
finite polynomial degree. This enrichment is intentional: ordinary Q2
pullback alone does not reproduce Cartesian quadratics on this torus.

Define the full conforming space and its original Dirichlet subspace by

    Vh = { v in C0(Omega closure;R^3) : v|K in U(K)^3,
           with identified periodic seam traces },
    Vh,0 = { v in Vh : v|_(r=a)=v|_(r=b)=0 }.

Continuity and zero trace mean equality of the **entire face functions**, not
merely agreement at the 27 Q2 nodes. This matters for the enrichment. A concrete
basis-independent assembly definition is: use the listed local generators,
quotient local linear dependencies, and impose equality of all coefficients
in the combined face trace spans on every shared or periodic face. The same
trace-span coefficient constraints impose zero on each radial face. These are
finite linear constraints defining a finite-dimensional element space. No
claim is made that the ordinary Q2 nodal degrees of freedom are unisolvent for
the enrichment. The ordinary globally conforming Q2 mapped space is a subspace,
so the constraints cannot remove its approximation functions.

Now specify the **Cartesian STF**, generally discontinuous tensor element. Let
E1=diag(1,-1,0), E2=diag(1,1,-2), and E3,E4,E5 have respectively the xy,xz,yz
symmetric pairs equal to 1 and all other entries zero. These span STF(3).
Let B_K=DF_K, g(qhat)=B_K^-T grad_xi qhat, and for Cartesian unit vector ei set

    D_i(qhat) = ei tensor g(qhat) + g(qhat) tensor ei
               - (2/3) (ei dot g(qhat)) I.

Define its reference tensor generators explicitly by

    That(K) = span{ D_i(xi1^alpha1 xi2^alpha2 xi3^alpha3) :
                    i=1,2,3; 0<=alpha1,alpha2,alpha3<=2 }
              + span{ E_beta p(F_K) : beta=1,...,5; p in P1(R^3) },
    T(K) = { That composed with F_K^-1 : That in That(K) },
    Th = { T in L2(Omega;STF(3)) : T|K in T(K) }.

This is a componentwise Cartesian composition map with the displayed
geometry-dependent gradient generators, **not** a Piola transform and not an
assertion that plain polynomial STF pullbacks contain the strain. The
polynomial inputs to the first generator have Q2 degree (their reference
partial derivatives have degrees (1,2,2) and permutations); the second has
Cartesian total degree 1. Inverse geometry factors and trigonometric functions
are part of the tensor mapping. This gives a finite, explicit space without
inventing a polynomial degree for nonpolynomial tensor components. There is no
interelement tensor trace matching: this auxiliary strain is L2 conforming,
not H(div) conforming. Requiring H(div) of this strain would be a different
specification.

## 3. Conformity and strain containment — PROVEN

Each local displacement function is smooth to its faces. Full matching traces
make the distributional first derivative equal the elementwise derivative:
interior face delta terms cancel. Hence Vh is a subspace of H1(Omega;R^3).
The periodic Cartesian traces agree even across a chart seam, because X and its
coordinate frame are periodic; no vector-component rotation is introduced.
Zero trace on both physical radial boundaries gives Vh,0 contained in
H1_0(Omega;R^3). There are no other physical boundary components.
Representation by actual functions, with redundant coefficients quotiented,
is faithful; coefficient redundancy cannot be mistaken for a displacement
kernel.

For the mapped Q2 part, the chain rule gives

    grad_x(qhat composed with F_K^-1) composed with F_K
        = B_K^-T grad_xi qhat.

Consequently L(ei qhat composed with F_K^-1) is exactly the generator D_i.
For the Cartesian P2 enrichment, its strain is a Cartesian STF polynomial of
total degree <=1 and lies in the second generator span. By linearity,

    L(Vh)|K subset T(K),   and therefore L(Vh) subset Th.

This is equality of weak strain with the represented tensor function, not
projection accuracy or a nodal product-rule assertion.

## 4. Exact integration, adjoint, elimination, and reactions — PROVEN

All pairings below are exact Lebesgue integrals on the exact geometry. In
particular, use the full Cartesian Frobenius product T:U=sum_ij T_ij U_ij:
off-diagonal entries are counted twice. In a five-component basis its actual
Gram matrix must be used. Set

    (T,U)_Omega = sum_K integral_cube (That:Uhat) |det B_K| dxi,
    (v,w)_Omega = sum_K integral_cube (vhat dot what) |det B_K| dxi,
    |det B_K| = (r q) det DA_K.

Every surface integral uses |cof(B_K) n_reference| and the physical oriented
normal, giving c q at a radial face in (theta,phi) coordinates. These pairings
are positive on nonzero represented functions. They are not conformal-volume
psi^6 pairings: the conformally rescaled momentum equation div LW=f uses
Euclidean dV; psi^6 belongs in f. The audit's physical volume weight is a
separate diagnostic.

**Integration requirement:** exact integrals define this proof-only method.
No finite tensor-product Gaussian order is claimed exact: inverse geometry
factors and trigonometric coefficients make the integrands nonpolynomial, and
psi-dependent loads are nonpolynomial too. Any future quadrature replacement
needs a new proof of positivity, uniform norm equivalence and controlled
variational error (including volume/face compatibility). Merely choosing a
rule exact on polynomial degree four does not satisfy this specification.
There is no implemented or certified finite quadrature algorithm here.

For symmetric trace-free T and v smooth enough for normal traces,

    (div T,v)_Omega = -(T,Lv)_Omega/2
                     + <T n,v>_boundary.

Indeed T_ij partial_j v_i = T:Lv/2 by symmetry and zero trace; the divergence
theorem gives the result with the stated weights and normals. Density gives
the same identity for v in H1 and T in H(div;STF), with the normal trace in
H^-1/2. The proof does not require an exact Leibniz rule on a finite nodal
algebra. For an arbitrary discontinuous Th tensor, apply the identity
cellwise instead: on an internal face with chosen normal n from K- to K+,
the boundary contribution is v dot (T^- - T^+)n. Periodic seam contributions
are paired in exactly the same manner. They do **not** generally cancel for
a discontinuous tensor. Thus the variational divergence includes the jump
contributions; it is not just the broken classical divergence.

For prescribed f in L2 and original homogeneous W data, define the mixed problem

    (sigma_h,tau_h) = (LW_h,tau_h)          for all tau_h in Th,
    -(sigma_h,Lv_h)/2 = (f,v_h)            for all v_h in Vh,0,
    (sigma_h,W_h) in Th x Vh,0.

Containment permits tau_h=sigma_h-LW_h, proving sigma_h=LW_h exactly.
The tensor mass form is invertible on functions. Elimination gives the
original Galerkin problem, with the correct negative sign for div L:

    a(W_h,v_h) = -(f,v_h),   a(u,v)=(Lu,Lv)/2.

In any independent bases with tensor mass Ht, strain coefficient map Lh,
and homogeneous injection I, the interior operator is exactly

    M A_D = -(Lh I)^T Ht (Lh I)/2.

Here M is the exact displacement mass if an L2 Riesz representation of the
interior variational divergence is desired. This is the weak adjoint identity;
there is no claim of a collocated strong divergence matrix or strong boundary
row moments.

Boundary reactions are retained, not set to zero. For a trace g_h from Vh,
choose any v_h in Vh with that trace and define

    R_h(g_h) = (f,v_h) + (sigma_h,Lv_h)/2.

Two such lifts differ by an element of Vh,0, so the discrete equation proves
lift independence. For the exact solution, the same functional is
<LW n,g> on both radial boundaries with their opposite normals. Interior
jump terms are accounted for by the volume variational functional; an
individual element's raw sigma_h n is not automatically the equilibrated
boundary reaction. This distinction avoids imposing an unauthorized zero
traction. Adjoint and discrete reaction identities are exact.

## 5. Dirichlet kernel, Korn bound, and auxiliary stability — PROVEN

For v in C_c^infinity(Omega;R^3), integration by parts yields

    integral partial_j v_i partial_i v_j = ||div v||_L2^2,
    ||Lv||_L2^2 = 2||grad v||_L2^2 + (2/3)||div v||_L2^2.

The second equality follows by expanding
|grad v + grad v^T -(2/3)I div v|^2. Density extends both equalities to
H1_0(Omega;R^3). This uses zero full vector trace on the entire physical
boundary, as in the original problem, not an assumption of zero normal
strain. Coordinate seams introduce no boundary terms on the physical domain.
By Poincare with the fixed domain constant CP,

    ||v_h||_H1^2 <= (1+CP^2)||grad v_h||_L2^2
                 <= (1+CP^2)||Lv_h||_L2^2/2.

The constant is independent of h, enrichment dimensions, and angular mode.
Thus ker(L|Vh,0)={0}; no nonzero represented old radial or checkerboard mode
can have zero strain. This is an H1 bound, stronger than a statement about a
particular coefficient norm. The variational momentum operator is uniformly
coercive after changing its sign, so the fixed-load Galerkin solution exists
uniquely. Also

    ||W_h||_H1 <= C ||f||_H^-1,   ||sigma_h||_L2 <= C ||f||_H^-1.

These follow from testing with W_h, Korn, and continuity of L. Hence auxiliary
elimination is uniformly stable in the physical H1 x L2 norms. No
mesh-independent Euclidean condition-number claim for a particular basis or
nonlinear coupled-system stability claim is made.

## 6. Reproduction, approximation, and boundary consistency

**Polynomial reproduction — PROVEN.** For every Cartesian vector polynomial
p of total degree <=2, p|K lies in U(K)^3 and Lp is represented exactly in
T(K). Such p also belongs to the unconstrained global Vh: it is single-valued
on all physical faces, including periodic seams. Thus constants, affine
vectors and Cartesian quadratics have exact strain reproduction, and the
pointwise weak derivative chain rule above commutes with their representation.
For any admissible test v_h,

    -(Lp,Lv_h)/2 = (div Lp,v_h) - <(Lp)n,v_h>.

This proves exact polynomial weak residuals with exact polynomial loads and
appropriate boundary lifts/reactions. These polynomial extensions are tests
of local consistency, not nonzero polynomials asserted to obey homogeneous
radial data. No strong collocated boundary formula is defined, so the former
C x^2=2 and cubic row tests cannot be attributed to this method.

**Approximation — PROVEN with stated regularity.** Let w be a periodic
Cartesian H3 field on Omega with zero radial trace. Pullback by X gives a
periodic H3 field on the fixed parameter shell. The conforming tensor-product
Q2 interpolant (or its bounded quasi-interpolation equivalent) respects the
matching periodic faces and homogeneous radial trace. Reference polynomial
approximation and scaling on the regular parameter mesh, followed by the
fixed smooth geometry norm equivalence, give a member I_h w in Vh,0 with

    ||w-I_h w||_H1 <= C h^2 ||w||_H3,
    ||w-I_h w||_L2 <= C h^3 ||w||_H3,
    ||Lw-L I_h w||_L2 <= C h^2 ||w||_H3.

For H3 functions nodal interpolation is legitimate in three dimensions; the
zero trace interpolates to the zero polynomial on a radial face. The enriched
space contains this Q2 interpolant. Constants depend on the fixed geometry and
mesh regularity, not h. Even without H3 regularity, density and the Q2 subspace
give H1 best-approximation tending to zero for H1_0 fields. Smooth STF tensors
also admit L2 approximation in Th because it contains elementwise Cartesian
P1 STF: local Taylor approximation on uniformly regular cells gives O(h^2)
for H2 tensors. No regularity of the actual nonlinear solution is inferred.

For the **fixed exact load**, Galerkin orthogonality and coercivity give

    ||W-W_h||_H1 <= C inf_(v_h in Vh,0) ||W-v_h||_H1,
    ||LW-sigma_h||_L2 <= C ||W-W_h||_H1.

Thus O(h^2) energy/strain convergence holds if the momentum solution is H3;
convergence without that rate follows from density. Perturbing the load adds
C||f-f_h||_H^-1. This is consistency and convergence for the original linear
Dirichlet momentum problem at fixed psi, not a proof for the coupled system.
No enhanced L2 rate for the computed W_h is assumed.

For discrete traces g_h with lifts v_h, exact and discrete reaction functionals
satisfy

    |R_h(g_h)-<LW n,g_h>|
       <= ||f_h-f||_L2 ||v_h||_L2
          + (1/2)||sigma_h-LW||_L2 ||Lv_h||_L2.

This perturbed-load statement assumes L2 loads. An interior H^-1 norm alone
does not control nonzero-trace boundary lifts. In particular, for the fixed
exact load the first term is zero. Taking
v_h=I_h v for each smooth periodic lifting v gives bounded lift norms and
convergence of the reaction against its interpolated trace to <LW n,trace v>.
This proves reaction consistency for these tests. It neither proves L2
convergence of raw face tractions nor a uniform reaction estimate on the whole
H1/2 trace space without a stable discrete lifting theorem for that norm.

**Boundary consistency — PARTIAL:** exact essential conditions, exact weak
boundary terms, convergence in H1 for the original linear Dirichlet problem,
and the reaction consistency just stated are established. Pointwise PDE
truncation at endpoints/transition layers, strong boundary-band residual
orders, and pointwise traction accuracy remain unresolved. A weak Galerkin
row is a functional on a test space, not an inaccurate strong stencil whose
boundary values may simply be declared consistent.

## 7. Unchanged Hamiltonian and frozen audit: exact obstruction and STOP

The finite-element scalar candidate needed to *state* the unchanged coupling
is psi_h in 1+Uh,0, where Uh,0 is the scalar conforming space from Section 2
with zero radial trace, restricted to strictly positive functions. For each
fixed mesh these functions and sigma_h are piecewise smooth and bounded, so
|sigma_h|^2 and all powers of a strictly positive psi_h are integrable. The
only permitted Hamiltonian strain contraction is the full Cartesian
|sigma_h|^2=|LW_h|^2; a five-entry unweighted sum, a projected square, or altered
physical E would change the problem. Exact integration would give the
original scalar weak equation with this contraction. This states the coupling;
it does not establish existence, positivity, uniform nonlinear estimates, or
convergence of a coupled solution.

There is some useful control, but it is insufficient. If sigma_h converges to
LW in L2 and is bounded there, then

    || |sigma_h|^2-|LW|^2 ||_L1
       <= (||sigma_h||_L2+||LW||_L2)||sigma_h-LW||_L2.

Consequently the quadratic source converges against bounded scalar tests. In
three dimensions L1 alone does not provide a uniform H^-1 bound or convergence
against every H1_0 test. For example, a smooth tensor bump with amplitude
O(epsilon^-1) and support radius epsilon has L2 norm O(epsilon^1/2) tending to
zero, while its squared magnitude has H^-1 norm of order epsilon^1/2, much
larger than its L1 norm O(epsilon). Increasing the tensor amplitude to
epsilon^-alpha with 5/4<alpha<3/2 still makes its L2 norm tend to zero, while
the square's H^-1 norm scales as epsilon^(5/2-2alpha), which diverges. This
scaling follows by testing against a compactly supported H1-normalized bump
of amplitude epsilon^-1/2. It illustrates the missing functional estimate,
not a counterexample claimed to be a solution of this discrete system.
Additional uniform integrability, for example convergence of sigma_h in
L^(12/5) giving squared-source convergence in L^(6/5), plus suitable control
and positivity of psi_h, would close this particular pairing issue. Those
estimates for the coupled problem have not been proved here.

The frozen [array audit](audit_run43_arrays.py) requires n=8,16,32,64,
psi.shape=(n+1,n,n), W.shape=(3,n+1,n,n), exact radial values, finite arrays,
and strictly positive psi. Its unchanged diagnostic
[boundary audit](../reproduce/diagnose_toroidal_boundary_convergence.py) has
Git blob `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22` at the starting commit.
The twelve checks are the two fixed boundary bands, the flux Hamiltonian and
second-order tensor momentum RMS metrics, and three refinement pairs; each
observed order must be finite and >=1.5. Fourth-order momentum and endpoint
quantities are additional diagnostics, not replacement acceptance criteria.

There is a specific faithful output map available in principle: evaluate the
continuous Cartesian W_h and psi_h at the audit's original toroidal nodes,
with no smoothing or modified radial values. Tensor values need not be chosen
at element interfaces because the audit recomputes its tensor from the W
array. This sampling map respects W=0 and psi=1 exactly and is single-valued
for continuous fields. It does **not** make the audit's reconstructed tensor
equal sigma_h. No different tensor, norm, band, or threshold may be substituted.

Here is the obstruction to the required audit proof. The momentum audit
composes two finite differences of sampled W; the Hamiltonian uses second
scalar differences and squared first differences. The established H1 error
and L2 strain error control neither sampled second differences nor their
boundary-band RMS. Indeed a Q2 function supported on a fixed number of adjacent
cells of diameter h, with amplitude A_h, has H1 norm of order
A_h h^(1/2) in three dimensions, but a typical nonzero sampled second
difference has magnitude of order A_h/h^2. Such functions exist already in
the unenriched conforming subspace, including adjacent to the radial boundary
with zero radial trace, taking the mesh and audit spacings comparable and
choosing a vertex basis function with nonzero nodal samples. There is no
uniform bounded map from its H1 norm to
these strong residuals. Even an O(h^2) energy error does not imply a decaying
strong-residual RMS with order >=1.5: take A_h=h^(3/2), giving H1 size
O(h^2) but second-difference RMS size O(h) on a fixed band containing the
patch, since its sampled volume is O(h^3). This argument
identifies a missing estimate; it does not assert that Galerkin errors must
realize these worst cases.

Compatibility would require, at minimum, a proof that the *actual coupled
solutions* are positive and sufficiently regular, a uniform nonlinear
source-control argument, and sampled-error/consistency estimates for precisely
the unchanged finite-difference compositions, including the one-sided radial
closures and fixed bands. Even asymptotic O(h^p) residual bounds alone do not
prove all three finite-resolution ratios at n=8,...,64: leading terms,
remainders and possible cancellation must also be controlled. No such result
follows from the spaces, Korn bound, or weak polynomial reproduction.

**Audit compatibility — UNRESOLVED. This required proof cannot be completed
from the analytical bounds and existing evidence used here. STOP at this
obstruction.** No alternative design, stronger unproved rate, numerical
experiment, solver construction, physical solve, GMRES, sweep, or Actions run
is introduced to bridge it. The Hamiltonian and audit are preserved, not
certified compatible merely because they can be written or sampled.

## Report and next step

| Obligation | Status |
| --- | --- |
| Spaces | SPECIFIED |
| Strain containment | PROVEN |
| Adjoint identity | PROVEN (variational, with exact integration and reactions) |
| Korn stability | PROVEN (fixed-domain H1 bound, independent of mesh) |
| Boundary consistency | PARTIAL |
| Audit compatibility | UNRESOLVED |

Only this specification and derivation are added. Repository validation for
this change consists of document review, diff/preservation checks, parent
verification, and remote branch SHA verification; no numerical code is run.
The pushed commit SHA is reported outside this document to avoid a
self-referential commit hash.

**One next necessary step:** separately authorize a proof-only analysis of
whether this same pair admits coupled regularity and sampled-error estimates
sufficient for the unchanged twelve-check audit, starting from the obstruction
in Section 7. No physical solves are authorized. STOP.
