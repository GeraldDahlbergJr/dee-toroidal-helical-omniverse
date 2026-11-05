# Run 43: one positive diagonal boundary norm feasibility study

Base: `fab11b7ed4defc34687d88ab7fdb277ee91998a7`. One agent; one
mathematical investigation. The latest six-row certificate, derivation and
report were read, not rerun. No historical artifact or physical source is
changed. This study stops at an exact inconsistent necessary subsystem.

## Fixed family and exact weighted SBP conditions

Use the same 16 dimensionless nodes x_i=i and d=h D_r. Rows 0,...,5
and their reflected outer rows may vary; rows 6,...,9 retain precisely
d_i,i-1=-1/2, d_i,i+1=1/2. All fifteen skew couplings within the six-node
boundary block remain allowed. No other couplings are introduced.
Relax only the radial norm to H=h W with
W=diag(w0,...,w5,1,1,1,1,w5,...,w0), every wj>0.
The physical weighted SBP equation is equivalently

    W d + d^T W = B,   B=diag(-1,0,...,0,1).

Thus Q=W d obeys Q+Q^T=B: Q00=-1/2, Q15,15=1/2,
other diagonal entries vanish, and off-diagonal entries are skew.
The fixed transition gives Q56=1/2, Q65=-1/2 and no other
boundary-to-interior coupling. The outer block is the negative reflection.
Consequently d00=-1/(2w0), d_ij=qij/wi and d_ji=-qij/wj
inside the left block, and d56=1/(2w5). These are all weighted
SBP conditions, including the transition. Centered adjacent rows force equal
interior weights. Even allowing their common value rho initially would give
rho=1: summing the left constant equations gives -1/2+rho/2=0.

Let A=Q[0:6,0:6], y=(0,1,2,3,4,5)^T, e5 the last unit vector,
b=e5/2, and w=(w0,...,w5)^T. The complete constant and affine
conditions for every variable row are

    A+A^T=-e0 e0^T,
    A 1+b=0,
    A y+6b=w.

The unchanged transition and centered rows differentiate constants and affine
polynomials exactly. Reflection supplies every outer equation. The new exact
certificate lists all twelve scalar constant/affine equations and all six
quadratic error numerators, without assuming any previously eliminated
trapezoidal parameterization.

## All affected boundary and transition quadratic moments

As in the latest certificate, the frozen normal principal adjoint composition
must satisfy d^2 1=0, d^2 x=0, d^2 x^2=2 at the unconstrained momentum
rows. The first two follow from d1=0 and dx=1. Set

    t=A y^2+36b-2 diag(w)y,
    eL_j=t_j/w_j (j=0,...,5), eL_j=0 otherwise.

Here y^2 denotes componentwise squaring. The necessary quadratic conditions
at every affected inner boundary/transition row are exactly

    R_i=(1/w_i) sum_{j=0}^5 A_ij t_j/w_j=0, i=1,...,5,
    R_6=-t_5/(2w5)=0.

Endpoint momentum row 0 is excluded by Dirichlet data. Rows 7 and 8 are
unaffected. Reflection gives rows 9,...,14, with the same zero conditions;
row 15 is excluded. Thus no transition equation has been dropped. These
constant, affine and quadratic equations are the complete necessary moment
system used in the previous test, now with the changed norm and coefficients.
Any higher accuracy demand includes this subsystem.

For a purely polynomial certificate put P=product_j w_j and define

    P_i=sum_j A_ij t_j product_{k!=j}w_k, i=1,...,5,
    P_6=-t_5 product_{k=0}^4 w_k/2.

Since weights are positive, P_i=0 is equivalent to R_i=0; specifically
P_i=P w_i R_i for i<=5, and P_6=P R_6.

## Exact contradiction

The constant condition and SBP give 1^T A=-e0^T+e5^T/2.
Multiplying the affine condition by 1^T and by y^T gives respectively

    sum_j w_j=5/2+3=11/2,
    sum_j j w_j=y^T A y+15=15,

because y^T A y=-(y_0)^2/2=0. Therefore

    sum_j t_j=25/2+18-2*15=1/2.

This is the same nonzero weighted quadratic error as before, now proved
without specifying any individual boundary weight. The affine condition and
SBP also give, column by column,

    sum_{i=1}^5 i A_ij-3 delta_j5=-w_j.

The term -3 is exactly the centered transition row 6 contribution. Hence

    sum_{i=1}^6 i W_ii R_i = -sum_j t_j = -1/2,
    sum_{i=1}^6 i P_i = -P/2.

All six residuals cannot vanish. For positive weights P>0, so the polynomial
equations require 0=-P/2<0. This is a precise algebraic contradiction,
valid for every derivative coefficient in this single family. It requires
neither nonlinear elimination nor numerical parameter fitting. The executable
certificate uses exact SymPy rational linear elimination only and verifies
the identities symbolically. Its largest matrix has six rows and six columns;
the enclosing stencil has 16 nodes. No angular tensor grid is assembled.

The same necessity for the continuum STF strain as in the latest derivation
applies: a locally quadratic Dirichlet field gives an order-one principal
normal momentum defect, while smooth geometry and lower-order corrections
vanish with h. The invertible normal STF principal symbol cannot remove it.

## Stop and scope of conclusions

Positive weighted norm satisfying the simultaneous requirements: INFEASIBLE.
Boundary moments: FAIL. A positive norm by itself is possible, but cannot
satisfy these derivative/momentum conditions. This is a failure of their
simultaneous feasibility, not a claim that positive diagonal matrices do not
exist. No feasible member exists to undergo the conditional strain-kernel or
adjoint checks. Their statuses for this requested construction are UNRESOLVED;
the historical family's adjoint result is not repurposed as a new check.

STOP at this contradiction. No kernel search, adjoint construction, alternative
family, second attempt, physical solve, nonlinear solve, GMRES, grid sweep or
GitHub Actions follows. No claim of uniform coercivity is made.

Next necessary step: separately authorize a change to the fixed
centered-transition restrictions before another feasibility investigation.

Reproduce only this new certificate with
`python research/run43_positive_boundary_norm/check_exact.py`.
