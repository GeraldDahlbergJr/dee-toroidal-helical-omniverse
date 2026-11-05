# Run 43: four-row boundary-block strain, exact rejection

Base: a6334b7c5757e56748de6b23cd55101e1f9df317.
One agent, one candidate. No previous experiment rerun. No production changes.

## Candidate and necessary moments (derived before implementation)

Replace only the first four radial derivative rows and their reflected outer
rows. Retain the original trapezoidal radial norm H_r=h diag(1/2,1,...,1,1/2),
original radial Dirichlet values W=0, angular periodicity, original continuum
momentum equation div LW=8 pi G psi^6 S and Cartesian STF strain definition.
This is a wider boundary-block alternative to the base's endpoint-only closure.
Require scalar SBP H_r D+D^T H_r=diag(-1,0,...,0,1), exact derivatives of
constants and affine functions, and the centered derivative beyond the block.
These are defining restrictions of this candidate, not necessary restrictions
on every possible closure. Dimensionless nodes are j, and D=d/h.

The six free skew entries in the first four rows of H_r D, followed by
constant/affine moment elimination, give exactly this one-parameter block:

    d[0,:5] = (-1, 2t, 2-4t, 2t-1, 0)
    d[1,:5] = (-t, 0, 3t-1, 1-2t, 0)
    d[2,:5] = (2t-1, 1-3t, 0, t, 0)
    d[3,:5] = (1/2-t, 2t-1, -t, 0, 1/2).

Rows j>=4 are centered until the reflected outer block. No coefficient fitting
or physical evidence is used. All coefficients and eliminations are rational.

For a frozen normal principal part the adjoint-induced interior operator is
D^2, by SBP and zero boundary trace. Necessary local polynomial conditions at
rows i=1,2 for vanishing truncation error are

    (d^2 1)_i=0, (d^2 j)_i=0, (d^2 j^2)_i=2.

The first two are automatic wherever the composition samples these rows.
The quadratic residuals are exactly

    i=1: 2(3t-2)(3t-1)=0,
    i=2: 4(2t-1)(3t-1)=0.

Their only common root is t=1/3. For second-order consistency the additional
cubic moments would be (d^2 j^3)_i=6i; their residuals are

    i=1: 6(2t-1)(5t-4),
    i=2: 6(2t-1)(11t-3).

At t=1/3 these are 14/3 and -4/3, respectively: second order already fails,
but these O(h) errors alone would still allow first-order consistency.
It is essential to check every affected row before claiming vanishing error.

## Exact contradiction: mandatory stop

The necessary first-two-row quadratic moments force t=1/3. At the next affected
interior row, however,

    (d^2 j^2)_3=37/18 != 2, residual=1/18.

At row 4, also affected through composition, the residual is -1/6.
Thus the necessary quadratic conditions on the boundary block are inconsistent:
first two rows force t=1/3, while row 3 requires 1/18=0. Investigation stops.
No closure implementation, second candidate, field solve, GMRES, nonlinear
iteration, grid sweep or physical comparison follows this contradiction.

This is also a genuine toroidal consistency obstruction. Choose a fixed smooth
radial scalar u equal to (r-a)^2 near the inner face, smoothly cut off before
the outer face, and W=u c for constant Cartesian c. Both original Dirichlet
conditions hold. Near the inner boundary, the leading principal radial term
is Q er u'', with Q er=c+er(er dot c)/3. Smooth J=r(1+r cos theta) and frame
coefficients change the quadratic principal-part error only by O(h_r), since
u=O(h_r^2) and D u=O(h_r) in the fixed-width block. Angular and lower-order
terms likewise contribute no O(1) cancellation. At theta=phi=0 and c=(1,0,0),
Q er=(4/3)c, so row 3 has limiting x momentum error 2/27. A nonzero limit
excludes every positive boundary consistency order. This argument does not
claim equivalence to the frozen audit or alter it.

## Adjoint, nonnegative energy, and coercivity limits

For this candidate's formal L, define the full weighted STF divergence by

    Div T = Hv^-1(-L^T Ht T/2+B T),
    B=E^T Hb N,  Hv=I_3 tensor H, Ht=I_9 tensor H,
    H=diag(J)(H_r tensor H_theta tensor H_phi).

This gives the exact original radial traction Green identity for every nodal
vector and STF tensor, with all nine tensor entries counted. No discrete
geometry product rule is presumed. For Dirichlet injection I, EI=0, hence
M A=- (L I)^T Ht (L I)/2, M=I^T Hv I. The strain energy is symmetric and
nonnegative for positive Ht. These are algebraic identities, even for this
inconsistent candidate; no original expanded discrete operator is replaced.

There is a further obstruction at the forced t=1/3: the first two interior
columns of d coincide. The radial array u=(0,1,-1,0,...,0) therefore has
D u=0, including the unchanged rows and separated reflected outer block.
For angular-constant Cartesian W=u c all angular derivatives vanish too,
so L W=0 exactly on the actual toroidal geometry. Thus strict strain coercivity
in any norm positive on this nonzero Dirichlet vector is disproven, not merely
unproved. Nonnegative energy does not imply a discrete Korn estimate.
No further coercivity experiments are performed after the moment failure.

## Certificate and preservation

`check_exact.py` records the exact rational moment elimination and contradiction
above only. It is a symbolic obstruction certificate, not an implementation of
a closure or a new investigation after the stop. It imports no repository
solver, assembles no physical grid and runs no physical tests. All changes are
new files under this directory; prior evidence, frozen audit, sources and solver
remain byte-for-byte unchanged from the verified base.

Next necessary step: obtain separate authorization before considering any
other boundary strain candidate. Physical convergence tests remain unauthorized.
