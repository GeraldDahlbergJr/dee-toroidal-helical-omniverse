# Run 43 — one boundary strain investigation

Base: `b075780b0c30296f0e7fbd6fb2230d503e31b40e`.
Branch: `research/run43-boundary-strain`. One agent; one candidate.

Read the base's `RUN43_ADJOINT_CANDIDATE_20261105.md`, its exact script and
certificate, the necessary-adjoint and quadrature-obstruction derivations, and
the historical closure reports to distinguish this construction from the
November 3–5 tests. Those tests are not rerun. This is a new endpoint SBP strain
construction, not a direct Hessian repair or retained-stencil quadrature change.

## Unchanged continuum target and exact old defect

The target is div LW = 8 pi G psi^6 S, with
LW_ij = partial_i W_j + partial_j W_i - (2/3) delta_ij partial_k W_k.
The Hamiltonian, original matter sources, G=.001, psi=1 and W=0 on
r=a=1/20, b=1/4, and angular periodicity remain unchanged. The recovered
continuum formulation has the provenance stated in the base derivation; this
report makes no new manuscript-verification claim.

Let H_r=h diag(1/2,1,...,1,1/2). For the old endpoint row
(-3/2,2,-1/2)/h, the first two interior columns of -H_r^-1 D_r^T H_r
produce, for F=J u', J=r(1+r cos theta),

    row 1: (-F_0 + F_2/2)/(h J_1),
    row 2: (F_0/4 - F_1/2 + F_3/2)/(h J_2).

Their constant-flux coefficient sums are -1/2 and +1/4. Specifically, the
endpoint entries D_01=2/h and D_02=-1/(2h), multiplied by endpoint weight
1/2, supply -F_0/h and +F_0/(4h). The interior centered rows cannot cancel
them. For u=(r-a)(b-r), the endpoint strain is already exact, so this is a
transpose/weight mismatch, not an endpoint strain approximation error.
The resulting leading terms are -(b-a) Q er/(2h) and
+(b-a) Q er/(4h), where Q er=c+er(er dot c)/3. Boundary traction cannot
repair these rows since Dirichlet injection I satisfies E I=0.

## The single candidate: minimal diagonal-norm SBP strain

Keep the centered interior radial derivative and trapezoidal norm. Require
H_r D + D^T H_r = diag(-1,0,...,0,1), with an endpoint row supported on
nodes 0,1,2. The (0,0), (0,1), (0,2) entries force respectively
D_00=-1/h, D_01=1/h, D_02=0. Reflection gives the backward outer row.
Thus the candidate is uniquely determined within this specified minimal class:

    D_s u_0=(u_1-u_0)/h,
    D_s u_i=(u_(i+1)-u_(i-1))/(2h), 1<=i<=N-2,
    D_s u_(N-1)=(u_(N-1)-u_(N-2))/h.

Use G_i=er_i D_s + et_i D_theta/r + ep_i D_phi/q and construct L_s with
the same Cartesian STF formula above. Only the two boundary strain rows
change. No normal derivative or traction is prescribed as extra data.
The endpoint derivative is first-order accurate; the consistency of its
adjoint composition must therefore be tested, not presumed from SBP.

## Full weighted adjoint and energy identity

Let H=diag(J)(H_r tensor H_theta tensor H_phi), with positive periodic angular
weights, Hv=I_3 tensor H and Ht=I_9 tensor H. All nine tensor components count,
so off-diagonal terms count twice. Let B map STF tensors to the original
outward radial traction pairing: V^T B T=(EV)^T Hb (T n),
n=-er,+er, Hb having surface weights r q h_theta h_phi. Define

    Div_s T = Hv^-1 (-L_s^T Ht T/2 + B T).

For every nodal V and STF T this gives exactly

    V^T Hv Div_s T = -(L_s V)^T Ht T/2 + (EV)^T Hb(T n).

This is a full variable-coefficient algebraic proof; scalar SBP alone does
not imply a discrete geometry product rule. With M=I^T Hv I and L0=L_s I,

    A_s=-M^-1 L0^T Ht L0/2,
    v^T M A_s w=-(L0 v)^T Ht L0 w/2.

Hence M A_s is symmetric, and the strain functional w^T L0^T Ht L0 w/4 is
nonnegative. -A_s is positive semidefinite. These properties survive the
consistency failure below, but do not make this an acceptable approximation.

## Exact boundary consistency counterexample — mandatory stop

Use only the smooth mathematical test W=u(r)c, u=(r-a)(b-r), c constant.
It satisfies both original radial Dirichlet values. Its centered derivatives
are exact at interior nodes. At the inner boundary D_s u=u'(a)-h; at the
outer boundary D_s u=u'(b)+h. With F=J u', the candidate radial transpose
contributions at the inner first and second interior rows are exactly

    R1=(F(a+2h)-F(a)+h J(a))/(2h J(a+h)),
    R2=(F(a+3h)-F(a+h))/(2h J(a+2h)).

The continuum radial contribution is F'/J. Write z=cos theta. Exact errors:

    R1-(F'/J)(a+h)
      = [a^2 z+a-4h^2 z]/[2(a+h)(1+(a+h)z)],
    R2-(F'/J)(a+2h)
      = -2h^2 z/[(a+2h)(1+(a+2h)z)].

The outer formulas replace a+h,a+2h by b-h,b-2h and the numerator a by b.
Thus the first-row error tends to +1/2 at BOTH faces; the second-row radial
error is O(h^2). The old O(1/h) constant-flux defect is removed, but the new
O(h) endpoint strain error is differentiated with coefficient 1/(2h).
For general smooth u the leading first-row defect is -u''/4 times Q er.

At theta=phi=0, c=(1,0,0), Q er=(4/3)c. For the exact toroidal shell the
first-row x momentum error therefore tends to +2/3 at both faces. At interior
rows the strain itself is exact for this test. Angular weighted-transpose
terms sample only that radial layer; hence endpoint strain errors cannot
enter them. Periodic centered differentiation of the smooth geometry-weighted
strain converges with O(h_theta^2+h_phi^2) error. In the continuum the remaining
geometry terms combine with F'/J to give div LW. They cannot cancel the
nonzero radial defect in the joint small-spacing limit. Thus the full first-row
truncation error is O(1), with nonzero limit; its boundary order is zero.
The second row on this test has O(h_r^2+h_theta^2+h_phi^2) error. This is a
counterexample, not a claim of improved accuracy for all admissible fields.
Multiplication by positive smooth psi^-10 cannot cure the first-row defect.

**Candidate FAILED. Adjoint PASS. No second closure is derived or tested.**

## Coercivity assessment and bounded certificate

Weighted symmetry and nonnegative energy are proved for the full candidate.
Strict positivity requires ker(L_s I)={0}; a uniform coercivity assertion
requires a mesh-independent discrete Korn bound in a stated norm. Neither is
established for this toroidal operator. Continuum Korn cannot simply be
transferred to collocated centered angular differences, whose even-grid
checkerboard derivatives vanish.

The exact seven-node normal reduction has positive scalar Gram leading minors
3/4, 3/8, 5/32, 15/256, 9/256. Its vector strain energy is a sum of positive
multiples (4/3 for the normal component and 1 for each transverse component)
of these scalar derivative energies. Thus that reduction is strictly positive;
it is not a toroidal or uniform Korn proof. Coercivity status: PARTIAL.
No further coercivity experiments are pursued after the consistency stop.

Reproduce with `python research/run43_boundary_strain/check_exact.py`.
SymPy 1.14.0; exact rational/symbolic arithmetic; largest matrix 7x7; no
assembled angular grid. The script verifies scalar SBP, the normal Dirichlet
Gram identity, its minors, and the exact toroidal radial errors and limits.
The full STF Green identity and angular consistency argument are proved above.
No production imports, physical fields, linear or nonlinear solves, GMRES,
full-grid sweeps, Actions, or historical experiments are executed.
Only this new research directory is added; production, original sources,
frozen audit, and historical evidence remain identical to the base.

Next necessary step: separately authorize one different boundary strain
construction with a pointwise-consistent adjoint. Physical convergence runs
remain unauthorized.
