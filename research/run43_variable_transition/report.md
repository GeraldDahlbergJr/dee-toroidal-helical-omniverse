# Run 43: one generalized boundary-to-interior transition family

Verified starting commit: `51bfca8feb707f0716ec3efc23a486389e80eb5d`.
One agent, no parallel work. Existing endpoint, four-row, six-row, positive-norm,
adjoint-candidate and quadrature-obstruction reports/certificates were read.
No completed certificate was rerun. This is a mathematical feasibility study,
not a production discretization or a physical solve.

## Unchanged physical problem

Retain the original Einstein–DEE continuum problem, including

    Delta psi + psi^-7 |LW|^2/8 + 2 pi G psi^5 E(psi) = 0,
    div LW = 8 pi G psi^6 S,
    (LW)_ij = partial_i W_j + partial_j W_i
              - (2/3) delta_ij partial_k W_k.

All original matter sources, Hamiltonian coupling, G=.001, toroidal shell
a=1/20, b=1/4, angular periodicity and radial psi=1, W=0 remain unchanged.
There is no additional traction, normal-derivative or source condition.
The frozen audit and historical evidence are unchanged. The original expanded
production matrix is not asserted to acquire an adjoint identity; only the
formal strain/adjoint construction below is assessed.

## One finite-support family and its polynomial requirements

Let x_i=i, d=h D_r, and retain precisely the original trapezoidal norm
H_r=h W, W=diag(1/2,1,...,1,1/2). Use any symbolic node count N>=18
to separate the two boundary neighborhoods. No N-by-N matrix or physical
grid is assembled: the certificate's largest matrix is 8-by-8.

The six boundary rows 0,...,5 and the two transition rows 6,7 are variable.
Permit **every** skew coupling of Q=W d within nodes 0,...,7, including all
boundary-to-transition couplings and the transition coupling Q67. Impose
Q00=-1/2, Qii=0 for i=1,...,7. The outgoing edge is Q78=1/2,
Q87=-1/2, as required by the first unchanged centered row 8. All other edges
from this block to nodes >=8 are zero. Rows 8,...,N-9 retain the nearest
neighbor centered derivative; the outer block is the negative reflection.
Thus the maximum variable index separation is seven; the only outgoing edge
has length one. All coefficients are dimensionless and independent of h.

This genuinely relaxes the previous fixed rows 6 and 7: they can sample any
of nodes 0,...,7. It retains a centered far field. The result below concerns
this explicitly bounded family, not unrestricted noncentered interiors or
nonlocal operators. No second family is tested.

Require exact differentiation of constants and affine polynomials at every
row. For vanishing local normal momentum error require the necessary
constant/affine/quadratic composition moments at **every** unconstrained
boundary or transition row:

    d1=0, dx=1;
    (d^2 1)_i=0, (d^2 x)_i=0, (d^2 x^2)_i=2.

The quadratic composition condition is necessary even for positive first-order
pointwise boundary consistency. No exact first derivative of x^2 is assumed
at the variable rows, and no higher polynomial order is imposed. Failure of
this subsystem already excludes higher-order demands.

## Exact equations for all affected rows

Set A=Q[0:8,0:8], y=(0,...,7)^T, H=diag(1/2,1,...,1),
b=e7/2, and 1 the eight-component constant vector. Weighted SBP and the
complete constant/affine equations are

    A+A^T=-e0 e0^T,
    A1+b=0,
    Ay+8b=H1.

These start with all 28 independent skew coefficients; exact linear
elimination leaves 15 free parameters. In particular the transition rows
are not forced back to their old centered form by assumption. The JSON gives
the full eliminated block, all 16 linear equations and all quadratic equations.
The unchanged rows satisfy the same constant/affine requirements identically.

Define the localized first-derivative quadratic error

    e=H^-1(A y^2+64b)-2y

on rows 0,...,7, and extend it by zero beyond that block. Componentwise
squaring is meant by y^2. Since d1=0 and dx=1, the composition constant and
affine moments hold identically. The complete remaining inner equations are

    R_i=(H^-1 A e)_i=0, i=1,...,7,
    R_8=-e_7/2=0.

Thus the new variable transition rows 6,7 AND the first unchanged row 8 are
included. Rows 9,...,N-10 are unaffected. Reflection gives the outer rows
N-9,...,N-2. Endpoint momentum rows 0 and N-1 are excluded under the original
Dirichlet data. Reflection reverses the sign of the first-derivative error
and preserves the composition residual, so it supplies the same zero
conditions. All equations are polynomial of degree at most two in the
eliminated skew parameters; the certificate lists all eight explicitly.

## Weighted adjoint identity and energy, before feasibility decision

Construct the formal strain L with this radial derivative, the unchanged
centered periodic angular derivatives, the original toroidal Cartesian frame,
and exactly the original STF formula. Let

    Hgeom=diag(J)(H_r tensor Htheta tensor Hphi),
    J=r(1+r cos(theta))>0,
    Hv=I3 tensor Hgeom, Ht=I9 tensor Hgeom.

Positive angular quadrature and the trapezoidal radial norm make these positive
definite. All nine tensor entries are counted. With E the radial trace,
Hb the original positive surface quadrature and Ntensor the original outward
traction map T -> T n (n=-er,+er), define the formal divergence on STF tensors

    Div T=Hv^-1[-L^T Ht T/2+E^T Hb Ntensor T].

Multiplying proves for arbitrary V,T

    V^T Hv Div T=-(LV)^T Ht T/2+(EV)^T Hb Ntensor T.

This is an exact full variable-coefficient identity. It does not assume a
discrete geometry product rule or identify the transpose with the original
expanded production operator. For the original zero-trace injection I,
EI=0, L0=LI and M=I^T Hv I>0 give

    M A_D=-L0^T Ht L0/2,
    -v^T M A_D v=||L0 v||_Ht^2/2 >= 0.

The corresponding strain functional is ||L0 v||_Ht^2/4. This proves weighted
adjoint identity PASS and positive **semidefinite** energy for every formal
member. The certificate independently verifies the arbitrary-matrix weighted
Gram algebra symbolically in a small abstract example; the full STF identity
is proved by the displayed multiplication, not inferred from that example.
Strict positivity, a kernel statement and a mesh-independent Korn estimate
do not follow from these identities.

## Exact contradiction and mandatory stop

From SBP and A1=-b,

    1^T A=-e0^T+b^T,
    1^T H e=49/2+32-2 sum_{i=1}^7 i=1/2.

From the affine equation and y0=0,

    y^T A-8b^T=-1^T H.

The first term represents rows 1,...,7; the second is exactly the row 8
transition contribution. Consequently every parameter choice satisfies

    sum_{i=1}^8 i W_ii R_i
       =(y^T A-8b^T)e
       =-1^T H e
       =-1/2.

Since W_ii=1 for i=1,...,8, simultaneous vanishing of all eight quadratic
moments would require **0=-1/2**. The executable certificate verifies this
identity over all fifteen remaining symbolic parameters with exact rational
arithmetic. No nonlinear solve, fitting or numerical search is used.

The obstruction is not caused solely by fixing rows 6 and 7: allowing both
to vary leaves the invariant. The remaining fixed centered far-field interface
and trapezoidal SBP/affine requirements force the weighted error sum.
Algebraically the same identity has 7 replaced by a symbolic block endpoint;
this explains the present failure, without testing another family or claiming
a no-go theorem for arbitrary noncentered interiors.

For the original continuum problem this is a necessary local test: choose a
smooth scalar u=(r-a)^2 near the inner face, cut off before the outer face,
and W=u c for constant Cartesian c. It respects the original Dirichlet data.
The frozen normal STF principal symbol is c+er(er dot c)/3, invertible for
nonzero c. The quadratic radial composition defect is order one; smooth
geometry and angular/lower-order corrections are order h locally and cannot
remove it. This is a consistency argument, not a modified physical source
or a physical solve.

Boundary moments: **INFEASIBLE**. Adjoint identity: **PASS** for the formal
family. Dirichlet strain kernel: **UNRESOLVED** because no moment-feasible
member exists and the required stop precedes its investigation. Coercivity:
**UNRESOLVED**; semidefinite energy alone establishes no strict or uniform
bound. No kernel experiment or alternative candidate follows the contradiction.

Only new research files are committed. No production, audit, source or historical
file is modified. No physical solve, GMRES, full-grid convergence study,
parameter fitting or GitHub Actions is run. Reproduce only this new certificate
with `python research/run43_variable_transition/check_exact.py`.

One next necessary step: separately authorize reconsidering the centered
far-field derivative restriction before any further feasibility study.
