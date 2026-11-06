# Run 43 — bounded analytical counterexample search

Base: `3422c928ea73cb514725b76a264b65754b7af2ec`.
Read only `research/RUN43_DECISIVE_INVARIANT_SET_20261106.md`.
One agent; analytical work only; unchanged formulation.

**Admissible counterexample: NOT ESTABLISHED. Result: UNRESOLVED.**

## New analytical finding: a coupled large-amplitude test

Consider a fixed admissible mesh and a nonnegative interior mapped Q2 cell
bubble η, extended by zero, normalized by max η=1. This is an admissible
function by the prescribed space definition. For every finite t≥0,
v_t=1+tη belongs to I_h(1,1+t), with the required entire radial trace 1.
This specifies a class of admissible inputs, not a certified counterexample
mesh or a violating input.

Write K_h for the original linear vector Galerkin solution operator acting
on the momentum source multiplier: all original constants, signs and vector
source S are retained in K_h, so W_h(v)=K_h(v^6). This uses the original
sixth-power momentum dependence underlying the designated report's M^6
estimate. Define Z_j=K_h(η^j), j=0,…,6, and τ_j=LZ_j. Linearity gives the
exact identities

    W_h(v_t)=Σ_{j=0}^6 binom(6,j)t^j Z_j,
    σ_h(v_t)=Σ_{j=0}^6 binom(6,j)t^j τ_j,
    F_h(v_t)=|Σ binom(6,j)t^j τ_j|²/[8(1+tη)^7]
             +a(1+tη)^5+b(1+tη).

These are identities in the unchanged coupled load family; the Z_j have
not been evaluated or replaced by freely chosen vector fields.

Let E={y:η(y)=0}. On this fixed mesh, as t→∞,

    t^(-12) F_h(v_t) → 1_E |τ_6|²/8   in L¹,
    t^(-12) R_h F_h(v_t)(x)
      → C_h(x):=(1/8)∫_E g_h(x,y)|τ_6(y)|² dy.

Proof: t^(-6)σ_h(v_t)→τ_6 in L² by the finite polynomial identity.
Its square therefore converges in L¹. The factor (1+tη)^(-7) is bounded
by 1 and converges pointwise to 1_E; dominated convergence applies to
|τ_6|², and the L¹ square-convergence controls the remaining difference.
The a and b terms divided by t^12 tend to zero in L¹, since 0≤η≤1
and a,b are bounded. The fixed-mesh continuous kernel slice g_h(x,·)
is bounded, proving the scalar limit. No mesh-uniform estimate is used.

Consequently, C_h(x)<0 for an actual mesh and its exact Z_6 would imply
a violation for every sufficiently large finite t. Conversely, comparison
for all finite upper endpoints M on that mesh would require C_h(x)≥0
for every such bubble and x. This necessary condition concerns the
momentum-generated strain outside the bubble support, not arbitrary
positive loads. A nonnegative C_h(x) does not decide the finite-t sign.
If an upper endpoint M is fixed in advance, the limit alone supplies no
violation within t≤M−1.

## Decision and limits

The bounded investigation stops at this exact reduction. No sign of C_h
has been certified on a specific admissible mesh; no Z_j coefficients,
finite violating t, or negative scalar evaluation have been established.
Thus the six required counterexample ingredients have not been supplied.
The limit is a conditional diagnostic, not evidence that a counterexample
exists. It neither proves nor refutes the m=1 lower comparison, even on
this bubble family. No further proof search is undertaken.

For m<1, the lower condition remains R_hF_h(v)≥m−1 for every input in
I_h(m,M). Neither that condition nor a positive invariant interval is
resolved. Even a future strict negative response at m=1 would refute
only the lower-endpoint-1 comparison on intervals containing that input;
it would not by itself exclude every interval with m<1, establish failure
of an upper bound, or establish nonexistence of a positive solution.

**One recommended next step, requiring new authorization:** certify the
sign of C_h(x) on one explicitly specified admissible mesh for one bubble,
retaining the exact original vector solve; pursue a finite-t certificate
only if that sign is strictly negative.

Preservation: one new research branch from the stated base and this report
only. No numerical experiment, physical solve, GMRES, sweep, Actions run,
source modification, or archive scan. Validation consists of the analytical
limit proof and Git parent/file/remote-head checks. The verified remote SHA
is reported externally to avoid self-reference. STOP; await authorization.
