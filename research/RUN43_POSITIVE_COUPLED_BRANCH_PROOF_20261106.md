# Run 43 — positive coupled branch: bounded proof investigation

Checkpoint: `ab84ef67b21c7ed9f7bcbb06bea6d2e6bb60d12f`.
Branch: `research/run43-positive-coupled-branch-proof`.
One agent; proof only. Only `RUN43_COUPLED_ESTIMATES_PROOF_20261106.md`
was read, with Section 5 as the starting point. Its definitions and accepted
linear mixed results are retained; no historical reports or source files were
consulted. No equations, physical sources, geometry, Dirichlet data, spaces,
or integration rules are changed. Results below are new implications of the
source structure stated there, not verification of additional source facts.

## 1. Source structure and an upper-barrier test

Use the fixed three-dimensional shell, scalar trace 1, homogeneous vector
Dirichlet data, and the original equations as recorded in the input report:

    -Delta u = |LW|^2 u^-7/8 + a(x)u^5 + b(x)u,
    div LW = 8 pi G u^6 S,
    a=2 pi G A >=0, b=2 pi G B >=0.

Assumptions for the following continuum calculations: G>0; bounded smooth
coefficients and S on the fixed smooth domain; admissible positive solutions;
and sufficient weak integrability for every displayed pairing. The accepted
Dirichlet Korn coercivity fixes the vector problem. No smallness of the physical
coefficients is assumed. The sign of the momentum weak load is immaterial to
the norm estimates.

A positive constant cannot be an upper supersolution wherever the scalar
source is strictly positive: its negative Laplacian is zero. This rules out
that particular barrier construction, not a bounded branch.

A nonconstant continuum barrier has a precise closure requirement. Let
z solve -Delta z=1 with zero trace and Z=||z||_infinity. For input
1<=v<=M, solve the vector equation and then the scalar Poisson equation with
all nonlinear factors evaluated at v and its vector solution. Smooth-domain
vector W^(2,p) regularity for some p>3 gives

    ||LW(v)||_infinity <= C_S M^6.

Here C_S includes the original source norm and fixed-domain regularity
constant; this is exclusively a continuum estimate. Put
K=C_S^2/8, a_*=||a||_infinity and b_*=||b||_infinity. The resulting nonnegative
Poisson load is bounded by

    K M^12 + a_* M^5 + b_* M.

Consequently the continuum interval [1,1+t z] is invariant under this frozen
Poisson map if, with M=1+tZ,

    t >= K(1+tZ)^12 + a_*(1+tZ)^5 + b_*(1+tZ).             (1)

This is a sufficient condition, not a necessary one. Increasing t indefinitely
does not enforce it. No quantitative constants or physical-source smallness
establishing (1) are supplied by the permitted report. Under (1), the stated
regularity and continuous dependence, the map is continuous and compact in
C(closed domain), so Schauder would yield a continuum fixed point. This
conditional construction supplies neither (1) nor a discrete invariant set.

## 2. Why the proposed strain exponent does not close an upper bound

The exponent 12/5 closes the scalar H^-1 load pairing under a positive lower
bound, as already established in the input. A Poisson L^infinity bound requires
more information. In dimension three, an L^q load estimate with q>3/2 is a
sufficient route. For the squared strain term this asks for strain in L^p
with p>3, in addition to the scalar lower bound. The endpoint p=3 needs an
additional endpoint estimate; it is not covered by this assertion.

The gap is real at the level of these norms. Choose a smooth nonnegative
bump eta supported in the unit ball, an interior point x0, and
f_r(x)=r^(-5/2) eta((x-x0)/r). Then ||f_r||_(6/5) is constant. The
zero-Dirichlet continuum Poisson solution has value at x0 of order
r^(-1/2), by the local three-dimensional Green kernel; its bounded regular
part contributes only O(r^(1/2)). Thus an L^(6/5) load bound cannot itself
give an L^infinity bound. These bumps are not asserted to be actual coupled
sources or solutions. They only exclude this norm implication. A structural
estimate for the actual coupled source could still bypass it.

## 3. Exact discrete stopping obstruction

For clarity define R_h f in the unchanged homogeneous scalar space by
(grad R_h f,grad q_h)=<f,q_h>. The frozen scalar map is

    T_h(v_h)=1+R_h[|LW_h(v_h)|^2 v_h^-7/8
                         +a v_h^5+b v_h].

The vector problem is uniquely solvable for each bounded positive input by
the accepted coercivity. In particular, testing that problem and using
Poincare/Korn gives ||LW_h(v_h)||_2 <= C M^6 for 1<=v_h<=M,
with C independent of h. This does not establish any stronger strain norm.

The first unclosed assertion in an invariant-set proof is that T_h maps a
whole-function positive interval into itself, uniformly in h. The source is
nonnegative, but no result in the permitted report proves R_h f>=0 even for
these particular coupled loads. Conforming enrichment does not make
truncations admissible. Neither the continuum torsion barrier nor continuum
comparison transfers through R_h by conformity alone.

For this route, one would need a proved discrete comparison/boundedness
property on the actual load family, or a uniform error estimate transferring
an independently established positive bounded continuum branch. A sufficient
package would include a mesh-independent positive invariant interval,
||LW_h(v_h)||_p <= C M^6 for some p>3, and a compatible uniform Poisson
L^q-to-L^infinity estimate (q=p/2>3/2), with quantitative closure of the
resulting upper-bound inequality. These are sufficient assumptions only;
none is adopted as a fact. Even an L^(12/5) strain bound would leave the
upper-bound gap in Section 2.

On each fixed mesh the frozen map is continuous on any closed interval
bounded away from zero, since it consists of linear solves and smooth
finite-dimensional nonlinear integrals. Such an interval is compact and
convex. Brouwer would apply if invariance were proved. Continuity alone
does not provide invariance. A continuation proof likewise needs estimates
preventing escape to zero or infinity and an appropriate degree or local
branch argument; no such estimates or starting branch are established.

**STOP:** the actual-source, mesh-uniform positive bounded invariant-set
estimate is missing. This is a failure to close the proof, not a proof that
the desired branch does not exist. No further audit or physical investigation
is undertaken.

## 4. Outcome

| Requirement | Status |
| --- | --- |
| Positive lower bound for actual discrete branch | UNRESOLVED |
| Mesh-independent scalar upper bound | UNRESOLVED |
| Uniform L^(12/5) strain integrability | UNRESOLVED |
| Coupled discrete branch existence/continuation | UNRESOLVED |

Exact remaining obstruction: no mesh-uniform control of T_h on a positive
bounded whole-function interval for the original coupled loads; the available
L2 strain bound, and even the proposed L^(12/5) source-pairing bound, do not
supply that control.

One next necessary step: prove or refute the positive bounded invariant-set
estimate for this exact Galerkin map and original source family. Until that
estimate is available, fixed-point existence and mesh-uniform continuation
remain unjustified.

Validation is limited to document review, the one-file diff, checkpoint-parent
and unchanged tracked-file verification, and remote SHA verification. The
verified commit SHA is reported externally to avoid self-reference.
