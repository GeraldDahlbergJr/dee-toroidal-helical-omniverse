# Run 43: necessary discrete adjoint and boundary-flux conditions

## Recovery and scope

Recovered `research/run43-product-rule-obstruction-20261105` at exactly
`4f8bcbe8e2ad2e33ce4ac100eb70b6d2c0b4ac73`. Read the November 5
`RUN43_PRODUCT_RULE_OBSTRUCTION_20261105.md`, the original solver and frozen audit,
and the tracked coordinated operator and boundary reviews. The checkpoint's
three product-rule SHA256 entries verify. No lost temporary report was sought.
This report derives necessary conditions for a discrete Green identity for the
original momentum equation; it does not assert that every convergent scheme must
have an exact Green identity. No physical equation, data, audit, or historical
artifact is changed. No numerical solve, grid sweep, Actions, or PRD submission.

## Original operator and weights

The conformal momentum equation is

    A W = div(L W) = f,       f_i = 8 pi G psi^6 S_i,
    (L W)_ij = partial_i W_j + partial_j W_i - (2/3) delta_ij partial_k W_k.

The original finite-grid solver in
`reproduce/solve_toroidal_helical_coupled_constraints.py` uses the expanded
blocks `Aexp_ij = delta_ij Lap_h + G_i G_j/3`, restricted to radial interior
unknowns. The frozen audit instead uses `Aten = Div_raw L_h`, with
`(Div_raw T)_i = sum_j G_j T_ij`. These finite matrices are distinct; neither
may be silently substituted for the other. Here

    G_i = M_(er_i) D_r + M_(et_i/r) D_theta + M_(ep_i/q) D_phi,
    q = 1 + r cos(theta),      J_volume = r q.

Use positive conformal Euclidean quadrature H approximating `r q dr dtheta dphi`.
Let Hv = I_3 tensor H and Ht = I_9 tensor H, where the tensor inner product sums
all nine Cartesian entries (so symmetric off-diagonal entries count twice).
Restrict tensor arguments to symmetric trace-free (STF) arrays. L_h is built
from the stated G_i without invoking any nodal Leibniz rule. The frozen physical
residual `psi^-10 Aten W - 8 pi G psi^-4 S` is psi^-10 times the conformal
residual. Thus the adjoint below concerns the conformal operator and conformal
weights, not the audit's physical RMS weight `psi^6 r q`. If expressing the
same pairing with that physical residual, its test weight must include psi^10.

## Necessary Green and traction identities

Let E_v extract radial boundary vectors, N extract traction `T n`, and H_b
approximate positive surface area `r q dtheta dphi` separately on each face.
The outward normal is `-er` at the inner face and `+er` at the outer face.
There is no angular boundary contribution because angles are periodic. Set

    B = E_v^T H_b N.

Continuum symmetry and zero trace imply
`sum_ij (partial_j V_i) T_ij = (L V):T/2`. Consequently the precise discrete
condition for the same integration-by-parts identity is

    Hv Div_h T + (1/2) L_h^T Ht T = B T                  (1)

for every STF T. This is equivalent to
`<V,Div_h T>_Hv = -<L_h V,T>_Ht/2 + <E_v V,N T>_Hb` for every V.
It determines Div_h on the STF space:

    Div_h = -(1/2) Hv^-1 L_h^T Ht + Hv^-1 B.             (2)

Equation (2) is a necessary algebraic characterization, not a proposed solver
or a justification for changing the frozen audit. For Div_raw, test (1) as a
matrix equality on an STF injection P:

    (Hv Div_raw + L_h^T Ht/2 - B) P = 0.                 (3)

It is not enough to check a few particular W or a single momentum sum. A stronger
componentwise sufficient condition is
`H G_j + G_j^T H = E^T H_b M_(n_j) E`. It is not asserted necessary separately
for every j: (3) is the relevant tensor condition. Variable coefficients,
quadrature, basis rotations and periodic angular terms all remain inside these
matrix identities. The November 5 obstruction precludes replacing them with an
ordinary nodal product-rule argument.

Let I inject vectors zero on both radial faces, and L0=L_h I. Necessary for the
original expanded interior operator to possess this same Dirichlet strain-energy
identity is

    I^T Hv Aexp I = -(1/2) L0^T Ht L0.                  (4)

For the tensor operator, (1) implies (4) with Aexp replaced by Div_h L_h because
E_v I=0. Therefore its weighted interior matrix is symmetric negative
semidefinite. Strict negativity further requires `ker(L0)={0}` (a discrete
Korn condition); the adjoint identity alone does not exclude spurious modes.
For the original Aexp, (4) is an additional independent compatibility test, not
a consequence of its continuum equivalence to div L. Any full-space original
operator claimed to equal Div_h L_h must also satisfy
`Hv Aexp + L_h^T Ht L_h/2 = B L_h`. No such equality for the actual toroidal
original operator is certified here.

## Integrated source and boundary reaction

For a constant Cartesian vector c, derivative matrices must give L_h c=0.
Testing (1) by c then gives the exact conformal flux balance

    c^T Hv Div_h T = (E_v c)^T H_b N T.                 (5)

For T=L_h W and a full nodal PDE, the right side carries the integral of f.
W=0 kills the boundary *work* term in the energy identity; it does not set N T
or the net traction to zero. Thus no zero-net-matter-momentum condition follows
from Dirichlet W=0. In an interior-only Dirichlet discretization, define the full
residual/reaction `rho = Div_h L_h W - f`, with I^T Hv rho=0. The exact statement
is instead

    (E_v c)^T H_b N L_h W = c^T Hv f + c^T Hv rho.        (6)

Endpoint PDE rows were replaced by boundary data, so the reaction term cannot
be discarded or presented as a failed source condition. Equivalently boundary
reactions/fluxes can be defined by the discrete weak residual. Consistency of
that flux with physical traction still needs proof. Neither zero traction nor
W'''=0 is an allowed additional boundary condition.

## Exact four-node certificate and its limits

`check_run43_adjoint_flux_exact.py` imports no repository solver. It uses only
SymPy rationals and one symbolic positive-weight calculation. There are four
spatial nodes, not a four-node-per-axis grid. It reproduces the original
three-point one-sided endpoint and centered interior D with h=1. The normal
reduction sets G_x=D, G_y=G_z=0 and uses constant Cartesian geometry; it is an
algebraic certificate, not a sampled toroidal field or a physical candidate.
The tensor and vector matrices have 36 and 12 slots respectively, but still
only four nodes.

On this reduction the original stencil cannot satisfy scalar endpoint-only
summation by parts for any positive diagonal H: its (0,2) entry is

    (H D + D^T H)_(0,2) = -h0/2,

whereas the boundary-only matrix has zero there. h0>0 makes equality impossible.
This obstruction also persists for the same one-sided/centered pattern on larger
normal lines with at least four nodes. It excludes a specific diagonal-norm,
endpoint-trace factorization; it does not exclude dense norms, another closure,
other trace definitions, or cancellation in the full toroidal tensor identity.

With H=diag(1/2,1,1,1/2), the certificate verifies that the raw normal tensor
composition fails (3), then verifies the exactly defined adjoint (2) satisfies
(1) on every STF tensor basis column. It verifies (4) for that adjoint-induced
operator and (5) for all constant Cartesian vectors. All six leading principal
minors of the normal-reduction Dirichlet strain Gram matrix are positive
(recorded in JSON), so this small reduction has no Dirichlet strain kernel.
This proves no discrete Korn bound for the toroidal system. The deliberate raw
stencil defect is the verified mathematical finding, not a failed certificate.

Reproduce in a fresh output directory:

    python research/check_run43_adjoint_flux_exact.py /tmp/run43-adjoint-UNIQUE

Verification completed with SymPy 1.14.0: exact identities PASS, all six minors
positive, checkpoint and prior SHA256 manifest PASS, frozen audit blob unchanged
`5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.
Only this derivation, certificate, result JSON, verification record and hash
manifest are added on the new research branch. Existing physical acceptance
results are unchanged.

## Next unresolved condition

Exhibit or rule out geometry-consistent positive quadratures and radial traction
traces satisfying (3) and the original expanded-operator condition (4) for the
actual toroidal coefficients, with a discrete Korn/coercivity bound and a
consistent boundary reaction flux. The four-node certificate does not establish
these whole-system conditions or convergence. No physical solve is justified by
this result alone.
