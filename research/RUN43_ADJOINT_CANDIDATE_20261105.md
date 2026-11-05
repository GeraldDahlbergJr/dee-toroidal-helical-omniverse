# Run 43: adjoint-induced candidate and boundary consistency obstruction

Base verified locally before mathematical work: `d206b2263bd0577c63026e1a3cc6ced5bb27a34d`.
Separate branch: `research/run43-adjoint-investigation`. One agent. Research only.

## Original equations and variational derivation

The equation provenance remains `docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md`
(the recovered action transcription, not an independently inspected manuscript),
`reproduce/solve_toroidal_helical_coupled_constraints.py`, and the base's
`research/RUN43_ADJOINT_FLUX_NECESSARY_20261105.md`.
The Einstein equation is G_mu,nu=8 pi G T_mu,nu for the original three-scalar
DEE action with kinetic matrix K_AB and potential U. The maximal conformally
flat constraints, with physical matter held fixed, are

    Delta psi + psi^-7 |LW|^2/8 + 2 pi G psi^5 E(psi) = 0,
    A W = div LW = f = 8 pi G psi^6 S,
    E = (Pi^T K Pi + psi^-4 sum_i D_i^T K D_i)/2 + U,
    S_i = -Pi^T K D_i.

Here G=.001; psi=1 and W=0 at a=1/20 and b=1/4; both angles are
periodic. L is the Cartesian trace-free symmetric gradient. No source,
physical equation, boundary condition, or acceptance threshold is changed.
The conformal momentum residual becomes the physical residual on multiplication
by psi^-10. This investigation concerns the conformal equation, not an
alternative physical weighting of it.

For prescribed f, the continuum functional is

    F(W) = integral |LW|^2/4 dV + integral f dot W dV.

Its first variation on zero-trace variations is integral (f-div LW) dot V.
Use the ORIGINAL nodal first derivatives and toroidal Cartesian coefficients

    G_i = er_i D_r + (et_i/r) D_theta + (ep_i/q) D_phi,
    q=1+r cos(theta), J=r q,
    (L_h W)_ij = G_i W_j + G_j W_i - 2 delta_ij sum_k G_k W_k/3.

Choose the natural positive diagonal quadrature H with radial trapezoid weights
(1/2,1,...,1,1/2) times h_r h_theta h_phi J. This is an explicit research
candidate choice, not a claim that the continuum uniquely chooses quadrature.
Let Hv=I_3 tensor H, Ht=I_9 tensor H, and I inject radial Dirichlet unknowns.
All nine tensor entries are counted; off-diagonal entries therefore count twice.
The discrete functional is (L_h Iw)^T Ht (L_h Iw)/4 + w^T I^T Hv f.
Writing M=I^T Hv I and L0=L_h I, its Euler equation is

    A_D w = -M^-1 L0^T Ht L0 w/2 = I^T Hv f mapped by M^-1.

Thus the adjoint-induced operator is DERIVED without modifying production.
Its negative is the Hessian of the strain functional in the M inner product.

## Exact integration by parts, symmetry, and coercivity limits

Let B=E^T Hb N be the original radial traction pairing, with outward normals
-er and +er and surface weights r q h_theta h_phi. On symmetric trace-free
(STF) tensors define

    Div_adj T = Hv^-1 (-L_h^T Ht T/2 + B T).

For every vector V and STF T, by matrix multiplication alone,

    V^T Hv Div_adj T = -(L_h V)^T Ht T/2 + (EV)^T Hb N T.

Because EI=0, I^T Hv Div_adj L_h I = -L0^T Ht L0/2 = M A_D.
This proves weighted symmetry and negative semidefiniteness for the full
variable-coefficient toroidal candidate, without a discrete Leibniz assumption.
Euclidean symmetry of M^-1 times the Gram matrix is not claimed. Coercivity
means positivity of -A_D: strict positivity requires ker L0={0}; a uniform
mesh coercivity claim requires a mesh-independent discrete Korn bound in a
specified norm. Neither full-toroidal property is proved here. The certificate
proves positivity only for a seven-node Cartesian normal reduction, using exact
positive leading principal minors. Symmetry and coercivity status: PARTIAL.
Boundary traction choices cannot change these interior rows since I^T B=0.

## Precise boundary obstruction: stop

The original radial derivative is centered in the interior and has endpoint row
(-3/2,2,-1/2)/h. Consider the inner face and the smooth Dirichlet test vector

    W(r,theta,phi)=u(r)c, u=(r-a)(b-r), c constant nonzero Cartesian.

This is a mathematical consistency test, not a physical solve or replacement
source. Its angular nodal derivatives vanish exactly and D_r u=u' exactly,
including endpoints. Thus T=L_h W equals the continuum strain at all nodes:

    T_ij = u'(r) [er_i c_j + er_j c_i - 2 delta_ij (er dot c)/3].

Put F(r)=J(r,theta) u'(r) at a fixed angle. The radial contribution in
Div_adj T factors into the smooth tensor contraction with er and, at the first
two interior rows, has exactly the scalar coefficients

    R_1 = [-F(a) + F(a+2h)/2] / [h J(a+h)],
    R_2 = [F(a)/4 - F(a+h)/2 + F(a+3h)/2] / [h J(a+2h)].

These follow from the COLUMNS of D_r in -H^-1 D_r^T H, not its rows.
The angular terms are bounded for this smooth field under angular refinement:
periodic centered differences of smooth geometry-weighted T have bounded
limits, and r,q stay away from zero. The continuum div LW is also bounded.
Consequently

    lim_(h->0) h (A_D W - A W)_(a+h) = -(b-a) Q er /2,
    lim_(h->0) h (A_D W - A W)_(a+2h) = +(b-a) Q er /4,
    Q er = c + er (er dot c)/3.

Q er cannot vanish for nonzero c. At theta=phi=0, c=(1,0,0), the x-component
limits are exactly -2/15 and +1/15 on the ORIGINAL shell. The certificate
computes these limits with the actual symbolic toroidal J=r(1+r cos theta).
This is a local Taylor argument for arbitrarily fine grids, evaluated exactly;
no large matrix or refinement sweep is built. No assumption of constant
Cartesian geometry is used in this boundary counterexample.

Hence the retained-stencil/trapezoidal-quadrature adjoint candidate is NOT
pointwise consistent with the continuum momentum operator at the first two
interior layers: the truncation error is O(h_r^-1), with nonzero leading term.
No positive boundary convergence order is established. Multiplication by a
smooth positive psi^-10 cannot remove the defect. The boundary values W=0 do
not require u'(a)=0; imposing that would change the authorized problem.
Away from the closure, centered radial and angular weighted-transpose terms
have the usual formal second-order consistency for smooth fields and geometry.
That bulk fact does not establish consistency on the whole Dirichlet space.
Continuum consistency status: PARTIAL (bulk formal consistency; strong boundary
consistency disproven). Boundary accuracy: O(h_r^-1) truncation defect.

This does not prove that the variational method cannot converge weakly in a
suitable norm: test functions vanish on the boundary and cancellations can
matter. Weak convergence and boundary reaction accuracy are unresolved. It also
does not exclude other quadratures, strain discretizations, or dense norms.
None is invented or tested after the obstruction. The required consistency
justification fails for this concrete candidate, so the investigation STOPS.

## Certificate and preservation

Reproduce only the bounded exact certificate:

    python research/check_run43_adjoint_candidate_exact.py

It imports SymPy and the standard library only; the largest matrix is 7x7.
It certifies a normal-reduction Green identity and strain Gram minors, then
computes the actual toroidal local boundary limits symbolically. The archived
JSON is `research/run43_adjoint_candidate_certificate.json`.
The full toroidal symmetry argument is algebraic above; no toroidal coercivity
claim is inferred from the normal matrix. Physical fields exist in historical
results, but the mandatory stop is reached before physical-field comparison.
No physical arrays were evaluated, no solve or GMRES was run, no full grid or
GitHub Actions was launched, and no previous experiment was repeated.
Only this new report, script, certificate, and their hash manifest are added.
Original production code, physical sources, frozen audit, historical evidence,
and thresholds are byte-identical to the verified base.

One next necessary step: separately authorize investigation of a boundary strain
and quadrature discretization with a consistent adjoint; production replacement
and physical solves remain unauthorized.
