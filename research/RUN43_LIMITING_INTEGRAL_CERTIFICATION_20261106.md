# Run 43 — limiting integral certification: inconclusive

Checkpoint: `e687be7c8be553f38ecfb51233ee4dabc549c66a`.
One agent; one attempt. The only research source read was
`research/RUN43_TARGETED_ANALYTICAL_COUNTEREXAMPLE_20261106.md`
at that checkpoint.

**Result: INCONCLUSIVE. No certified sign and no actual m=1 counterexample.**

## Specification obstruction

The permitted source gives
[
C_h(x)=\frac18\int_{\{\eta=0\}}g_h(x,y)|LK_h(\eta^6)(y)|^2\,dy,
]
but does not specify the physical domain, its coordinate map, mesh
admissibility constraints, full mixed finite-element spaces, scalar or vector
bilinear forms, boundary conditions for those operators, or the original
constants and vector source S. It explicitly leaves K_h as the original
operator with those data retained. Neither K_h nor g_h can therefore be
constructed from this source alone.

Exact mesh: **not specified or certified**. Bubble on a physical mesh:
**not specified or certified**. Physical evaluation point: **not specified
or certified**. These are unmet task requirements, not implicit choices.
Assigning arbitrary geometry, operators, or spaces to fill these gaps would
violate the preservation requirement. No mesh was tried.

For precision, on a known mapped hexahedral cell T=F([-1,1]^3), the usual
reference bubble would be
[
\widehat\eta(\xi)=\prod_{i=1}^3(1-\xi_i^2),\qquad
\eta|_T=\widehat\eta\circ F^{-1},\quad\eta=0\text{ off }T.
]
It is nonnegative, has maximum one, vanishes on every cell face, and belongs
to the mapped tensor Q2 space. This algebraic observation supplies neither
the missing physical cell map nor certification of membership in the
unspecified original mixed space. It is not an admissible-mesh certificate.

## What follows rigorously from the supplied reduction

Writing E={eta=0} and tau_6=LK_h(eta^6), boundedness of the fixed-mesh kernel
slice gives the conditional estimate
[
|C_h(x)|\le\frac18\|g_h(x,\cdot)\|_{L^\infty(E)}
                       \|\tau_6\|_{L^2(E)}^2.
]
This follows directly from the triangle inequality. The source supplies no
values or signed bounds for these factors; the estimate certifies no sign.
Nonnegativity of |tau_6|^2 alone cannot determine the signed kernel integral.
No numerical approximation or physical solve was performed.

If a future fully specified instance establishes C_h(x)=-delta<0, the
source's proved limit implies that for some finite T and every t>T,
[
t^{-12}R_hF_h(1+t\eta)(x)<-\delta/2.
]
For any such finite t, the source asserts 1+t eta belongs to I_h(1,1+t),
so the response is strictly negative on that admissible input. This is a
conditional finite-input implication, not a counterexample established here;
no negative delta, threshold, or finite violating input has been certified.
A nonnegative limiting integral would likewise not settle finite-t behavior.

## Stop and different strategy

Stop after this attempt. Recommend a specification-first algebraic audit:
supply a self-contained exact definition of the frozen geometry, admissible
mesh, mixed spaces, both discrete operators, and original source, then
analyze the resulting finite-dimensional positivity problem directly.
Do not commission another limiting-integral report with the same missing
operator data. This recommendation requires no change to the frozen model.

Preservation: only this new report on one new research branch from the
checkpoint. No equations, sources, spaces, geometry, frozen audit, or
existing files changed. No nonlinear solver, GMRES, sweeps, or Actions used.
Publication validation: check the commit parent, sole changed path, and
remote branch SHA; report the verified SHA externally.
