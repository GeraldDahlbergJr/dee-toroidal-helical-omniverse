# Run 43 — decisive invariant-set investigation

Checkpoint: `3e2b5fa50c68e289647594c1f46407aaf4994b91`.
Branch: `research/run43-decisive-invariant-set-20261106`.
One agent; proof only; exact geometry, integration, spaces and physical parameters.

| Requirement | Status |
| --- | --- |
| Discrete positivity on the actual coupled load family | UNRESOLVED |
| Mesh-independent upper bound for the frozen map | UNRESOLVED |
| Positive bounded whole-function invariant interval | UNRESOLVED |

No mathematical counterexample within the prescribed space **and actual load
family** has been established. The results below identify exact comparison
conditions and quantitative closure requirements, not nonexistence of a branch.

## 1. Scope and unchanged definitions

Read the requested positive-branch report. Consulted only its two directly
referenced reports for the space/source definitions and the three original
source/formulation files named there; no historical archive scan or prior-proof
reproduction. Root and research AGENTS.md were absent at the checkpoint.
The scalar space U_h,0 has local pullback Q2 + span{1,F_i,F_i F_j}, quotienting
dependencies and imposing equality of entire face traces, with zero radial
trace. The vector space is its Cartesian triple. Meshes have the specified
bounded aspect ratios and comparable sides within each brick; global
quasi-uniformity is not added. All statements concern represented functions.

Write I_h(m,M)={v in 1+U_h,0 : m<=v(x)<=M throughout the closed domain},
where 0<m<=1<=M. Let sigma_h(v)=LW_h(v), with the accepted variational
momentum equation, and

    F_h(v)=|sigma_h(v)|^2 v^-7/8 + a v^5 + b v,
    T_h(v)=1+R_h F_h(v),   (grad R_h f,grad q)=(f,q).

Here "/8" divides the strain term by eight; the inverse power is v^(-7).
No projection of the square or quadrature approximation is used.

The original source simplifies exactly. With d=(.25,.20),
P=(.12,-.04), K=[[1.2835,.0945],[.0945,1.189]], and
D=e_theta/r-3 e_phi/q, q=1+r cos(theta),

    A=(P^T K P)/2+(.945-1)^2 = .0127638,
    B=(d^T K d)|D|^2/2 = .068614375 (r^-2+9q^-2),
    S=-(P^T K d)D = -.030316 D,
    a=2 pi (.001) A,   b=2 pi (.001) B.

These decimal coefficients are exact algebraic expressions in the stated
decimal parameters. In particular |D|^2<=416, B<=28.54358,
|S|<=.030316 sqrt(416). Volume V=3 pi^2/25. Denote a_*=a and
b_*=2 pi (.001)(28.54358). These are unchanged, bounded, nonnegative sources;
inputs need not inherit the axisymmetry of the coefficients.

## 2. Exact discrete comparison criterion

For any independent real basis phi of U_h,0, put
H_ij=(grad phi_i,grad phi_j), which is positive definite. Define

    g_h(x,y)=phi(x)^T H^-1 phi(y).
    R_h f(x)=integral g_h(x,y) f(y) dy.                   (2)

This follows directly by solving the variational coefficient equations and
holds for every integrable load on each fixed mesh. It is basis independent.
For arbitrary nonnegative integrable loads, R_h f>=0 everywhere is equivalent
to g_h(x,y)>=0 everywhere: sufficiency follows from (2); necessity follows by
a nonnegative smooth bump supported where a negative continuous kernel slice
stays negative. This is an exact criterion, not an assertion of its sign.

For the actual family the necessary and sufficient criterion is instead

    integral g_h(x,y) F_h(v)(y) dy >= 0
      for every x and v in I_h(m,M).                    (3)

For interval invariance the exact lower condition is >=m-1, and the upper
condition is <=M-1. Thus failure of (3) would refute the lower endpoint 1
but need not refute every positive interval with m<1.

Split g_h=g_h^+-g_h^-. The actual comparison condition (3) is precisely

    integral g_h^+ F_h(v) >= integral g_h^- F_h(v).       (4)

The source has a positive floor a m^5+b m, but that floor alone supplies
neither side of (4) in the needed direction. A sign-changing kernel or a
nonpositive coefficient of H^-1 alone supplies no counterexample to (3):
the loads are tied to v and the vector solve. Conversely, source positivity
does not prove (3). Truncations are not admissible by conformity alone, and
no discrete maximum principle has been assumed.

No sign or actual-family compensation estimate for (4) is established here.
This is the first missing estimate on the route with lower endpoint 1.
For other positive endpoints the corresponding missing estimate is the
uniform lower bound m-1 in (2).

## 3. What the energy estimate actually gives

An explicit safe Poincare constant for this shell is

    C_P <= (.20/pi) sqrt(25/3).

Indeed J=rq lies between .0375 and .3125; apply the one-dimensional radial
Dirichlet Poincare inequality before these weight bounds, and use
|partial_r w|<=|grad w|. The same constant works componentwise.
Testing momentum and using ||grad W||_2<=||LW||_2/sqrt(2) yields

    ||sigma_h(v)||_2 <= C_2 M^6,
    C_2=8 pi G sqrt(2) C_P ||S||_2
       <=8 pi G sqrt(2) C_P (.030316 sqrt(416 V)).       (5)

Thus a genuine uniform actual-source estimate is available:

    ||F_h(v)||_1 <= m^-7 C_2^2 M^12/8
                    + V(a_* M^5+b_* M).               (6)

It does not imply uniform scalar supremum control.

For clarity define the exact finite-mesh constants

    D_h=max_x sqrt(g_h(x,x)),
    B_h=max_(x,y) |g_h(x,y)| = D_h^2.

The equality follows from Cauchy-Schwarz in the H^-1 coefficient inner
product and taking a diagonal maximum. Also

    D_h=sup_(q_h!=0) ||q_h||_infinity/||grad q_h||_2,
    ||R_h f||_infinity<=B_h ||f||_1.                    (7)

These are valid inverse/evaluation estimates for this enriched space,
without an unproved reference-basis conditioning claim. They are not uniform:
an interior mapped Q2 cell bubble, extended by zero, is already in U_h,0.
On a regular cell of size ell its amplitude-one value is 1 after normalization
and its gradient norm is <=C ell^(1/2). Hence D_h>=c ell^(-1/2) and
B_h>=c^2 ell^(-1). Refinement therefore prevents (6)-(7) from giving a
mesh-independent upper bound. This proves failure of this L1 operator-norm
route, not divergence for the actual coupled loads. No arbitrary bump is
claimed to equal F_h(v).

## 4. Quantitative upper closure and nonlinear growth

A sufficient alternative package would be, for some p>3 and q=p/2,

    ||sigma_h(v)||_p <= C_p M^6,
    ||R_h f||_infinity <= C_R ||f||_q,                  (8)

uniformly over the prescribed meshes and inputs. These are hypotheses,
not consequences of (5). Under (8),

    ||T_h(v)-1||_infinity
      <= alpha M^12 + beta M^5 + gamma M,
    alpha=C_R m^-7 C_p^2/8,
    beta=C_R V^(1/q) a_*,
    gamma=C_R V^(1/q) b_*.

If the lower comparison is separately available, a sufficient upper condition is

    alpha M^12 + beta M^5 + gamma M <= M-1.             (9)

A weaker family-specific scalar estimate could replace (8); neither is
necessary for existence. The momentum square produces the M^12 factor.
The common dependence of sigma_h and v might improve this overestimate,
but no cancellation or weighted estimate doing so is proved.

For positive beta, (9) requires M>1. Dividing by M shows that closure is
equivalent to finding M>1 with

    alpha M^11 + beta M^4 + gamma <= 1-1/M.             (10)

In particular necessary conditions for this sufficient test are

    gamma<1,   alpha<=11^11/12^12,   beta<=4^4/5^5,

by separately maximizing (M-1)/M^12 and (M-1)/M^5.
They are not jointly sufficient. An explicit sufficient test, for any chosen
M0>1, is alpha M0^12+beta M0^5+gamma M0<=M0-1.
Large M cannot rescue (9); its higher powers eventually dominate.
The physical a_*, b_* and C_2 are explicit above, but C_p and C_R in (8)
have not been bounded for the exact enriched spaces. Consequently (9) is
not certified for the unchanged physical parameters. Nor is its failure
established, and failure would not refute the branch.

An L^(12/5) strain bound would give only an L^(6/5) squared source.
It supplies neither estimate in (8) nor supremum control. No such inference,
continuum-to-discrete comparison transfer, or global quasi-uniform inverse
estimate is used.

## 5. Decision and preservation

The exact obstruction is the unestimated negative-kernel contribution in
(4) for the restricted momentum-generated loads, followed by a missing
mesh-uniform scalar bound for their weighted strain squares and quantitative
upper closure. The inverse evaluation calculation proves why the available
energy/L1 route cannot supply that bound. It is not a counterexample to
the actual map or branch.

**One necessary next step:** establish or refute the uniform lower comparison
in (2) for F_h(I_h(m,M)), for some fixed 0<m<=1<=M, starting with (4) for
m=1. Any refutation of an invariant interval must specify an admissible
whole-function input, its exact vector Galerkin solution, and the violating
scalar evaluation; arbitrary loads do not suffice.

Only this report is added. Validation: source algebra and analytical
derivations reviewed; remote commit parent, changed-file list, report contents
and branch-head SHA checked. No implementation, physical solve, GMRES, grid
sweep, fitting, regularization, source change or Actions run was performed.
The verified remote SHA is reported externally to avoid self-reference.
STOP; no additional investigation is authorized.
