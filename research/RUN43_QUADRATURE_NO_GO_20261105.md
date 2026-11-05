# Run 43: incompatible original expanded energy and nodal quadrature

Base: verified commit `3f294581597a1831edf6001d8a191dae134b1499`.
This is a new necessary-condition test only. No previous certificate, physical
solve, sweep, or Actions run is repeated. Only this report, its exact certificate,
and its output are added on a separate research branch.

## Decision and scope

**Compatible construction: NO**, retaining the original collocated gradients,
expanded momentum matrix, and positive nodal toroidal quadrature. In fact the
quadratic energy identity alone is impossible. This excludes arbitrary positive
nodal weights, including nonseparable geometry-dependent weights, and any
boundary traction consistent with zero boundary work for zero radial trace.
It does not exclude redesigning the discrete derivatives/operator or using a
nonlocal matrix inner product instead of a nodal quadrature. Such redesigns
would be a different task. No claim of impossibility for the continuum problem
is made.

The requested energy condition is at least

    W^T Hv Aexp W = -(L W)^T Ht (L W)/2,   E W=0,
    Aexp_ij = delta_ij Lap + G_i G_j/3.

The argument below tests only this quadratic equality; it does not assume the
stronger unsymmetrized matrix equality in the preceding report.

## Exact counterexample

Take five radial nodes on the original shell, with radii
`1/20, 1/10, 3/20, 1/5, 1/4`, and four nodes in each periodic angle.
Set `u=(0,0,1,0,0)^T` and `chi=(1,-1,1,-1)^T` around phi. For any fixed nonzero
Cartesian vector c, consider the two vector arrays

    W_plus(r_a,theta_t,phi_p) = u_a c,
    W_minus(r_a,theta_t,phi_p) = u_a chi_p c.

Both vanish on the original radial faces. With the original centered periodic
first derivative, `D_phi 1=D_phi chi=0`, whereas the direct second derivative
satisfies `D_phiphi 1=0` and `D_phiphi chi=-4 chi/h_phi^2`.
These statements follow from the rational 4-by-4 shift matrices in the
certificate; `h_phi=pi/2` remains an exact symbolic scale. The dimensionless
5-by-5 radial matrix has `D_r u=(-1/2,1/2,0,-1/2,1/2)^T/h_r`.
In particular `D_r u=0` on the only layer where W is nonzero.

At every node, for every Cartesian derivative and component,

    G_i W_minus,j = chi_p G_i W_plus,j.

The theta derivative of either scalar amplitude is zero; the phi first
derivative is zero; and the remaining radial term has precisely this sign.
Consequently `L W_minus = chi_p L W_plus` pointwise, and their Frobenius strain
norms, hence their weighted strain energies, are identical for every nodal
quadrature. This includes all boundary strain samples.

Only the middle radial layer contributes to either expanded work pairing.
On that layer `G_j W_plus,j=G_j W_minus,j=0`. In a second Cartesian derivative,
angular terms differentiate an array identically zero on that layer, so vanish.
The radial terms in `G_i G_j W_minus` equal chi times those in
`G_i G_j W_plus`; their work contributions therefore agree. This reasoning
applies componentwise and includes the grad-div cross terms, without a product
rule assumption. The radial Laplacian contributions also agree. The only
difference is the direct phi second derivative in the expanded Laplacian:

    E_minus - E_plus
      = -4 |c|^2/h_phi^2 sum_(t,p) H_(middle,t,p)/q_(middle,t)^2 < 0.

Here `q=1+r cos(theta)>0`. On the four theta nodes at r=3/20,
`q=(23/20,1,17/20,1)` is rational. Every H weight in the sum is positive.
The right-hand sides of the proposed strain-energy identities agree exactly,
while the left-hand sides differ strictly. At least one of these two explicit
Dirichlet vectors therefore violates that identity for every choice of positive
nodal quadrature. This is a joint counterexample family, not a choice of weights
fitted to a single failed vector. The same argument works on any even phi ring
with an interior radial spike whose centered derivative vanishes at its support.

## Consequences for traction and coercivity

Any tensor-adjoint construction satisfying the required Green identity has
boundary work `W_boundary dot traction=0` for both vectors. Changing its traction
cannot change the contradiction. Defining divergence by a weighted transpose
can satisfy the tensor-adjoint identity, but its induced operator cannot also
have the original expanded quadratic energy identity. Thus all four requested
properties cannot hold simultaneously; no candidate construction is advanced.

**Coercivity: UNRESOLVED.** The incompatible conjunction has no construction
whose coercivity can be proved. This counterexample does not show that the
expanded operator itself is noncoercive, nor disprove coercivity of a redesigned
adjoint operator. No such separate investigation is performed.

## Reproduction and preservation

Run `python research/check_run43_quadrature_no_go_exact.py`. The certificate
imports only SymPy and standard-library modules, uses rational matrices of at
most 5-by-5 and positive symbolic weights, and asserts strict negativity of
every weight coefficient in the exact work difference. It assembles no toroidal
grid matrix and solves no physical or algebraic field problem. The accompanying
JSON records the exact checks and formula. The Cartesian/toroidal identities
connecting these small matrices are proved above; no numerical pi evaluation
or floating-point tolerance is used.

Original sources, frozen audit, prior certificates, historical results, and
production branch are unchanged. The next necessary step, requiring new
authorization, is to decide whether the original expanded discrete operator may
be replaced by an adjoint-induced operator; quadrature and traction alone cannot
repair this obstruction.
