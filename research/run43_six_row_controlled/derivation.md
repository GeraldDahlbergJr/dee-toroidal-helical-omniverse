# One six-row controlled discretization family

Starting commit: 80c48400b80c8fd391f10cdd06a4f1d9c53ec490.
One agent. The latest feasibility report and its four-row, endpoint,
quadrature-no-go and necessary-adjoint exact certificates were read, not rerun.
The new certificate imports no production code. All matrices have at most
16 rows and columns; there is no assembled angular grid.

## Definition and complete necessary moment system

Set dimensionless nodes x_j=j, j=0,...,15, H=diag(1/2,1,...,1,1/2),
and physical derivative D_r=d/h. Only rows 0,...,5 and their reflected outer
rows vary. Rows 6,...,9 remain centered. Permit every skew coupling within
the first six nodes; the coupling d_56=1/2 is fixed by the centered row 6.
No other coupling is allowed. Require Hd+d^T H=B, B=diag(-1,0,...,0,1),
d1=0 and dx=1. This is a genuinely wider block than the rejected four-row
family, with six independent parameters after linear elimination. The script
starts with all fifteen independent skew entries, rather than assuming the
eliminated form. The certificate lists every resulting row coefficient.
The outer block is d_(15-i,15-j)=-d_ij.

Retain the continuum Cartesian STF strain, toroidal frame and positive
geometry weight J=r(1+r cos(theta)), original momentum equation
 div LW=8 pi G psi^6 S, Hamiltonian coupling, sources, G=.001,
a=1/20, b=1/4, psi=1 and W=0 on radial faces, and angular periodicity.
Use the existing centered periodic angular derivatives in the proposed strain.
This defines a mathematical family only, not an implemented physical operator.

For the frozen normal principal part, SBP makes the induced interior
adjoint composition d^2. The necessary moments for vanishing local error are
(d^2 1)_i=0, (d^2 x)_i=0, (d^2 x^2)_i=2.
The first two hold identically. Put e=d x^2-2x; split e into its inner and
outer supported parts. The inner part eL is supported in rows 0,...,5.
Every affected inner or transition row is i=1,...,6, with residual
 R_i=(d eL)_i = sum_(j=0)^5 d_ij [sum_(k=0)^6 d_jk k^2-2j].
This formula, together with the complete coefficient table and six expanded
polynomials in certificate.json, specifies all necessary quadratic conditions
R_1=...=R_6=0. Rows 7,8 are unaffected; rows 9,...,14 have the reflected
conditions. Endpoint momentum rows are not imposed under Dirichlet data.
Any higher requested polynomial order must include this inconsistent subsystem.

## Exact inconsistency, including the transition row

The script verifies identically over all six free parameters that
 1^T H eL=1/2,
 sum_(i=1)^6 i H_ii R_i=-1/2.
There is also a short SBP explanation of the second identity:
x^T H d eL=x^T B eL-(d x)^T H eL=-1/2.
Here x_0=0, the outer trace of eL vanishes, and dx=1. Only rows 1,...,6
contribute on the left. The first identity follows by summing the derivative
of x^2 up to the centered transition: its interface flux is
[m^2+(m+1)^2]/2 at m=5, while sum_(i=1)^m 2i=m(m+1);
the difference is 1/2. The exact script independently checks both identities.
Thus all six quadratic residuals cannot vanish, even over complex parameters.
No nonlinear equation solving, numerical fitting or previous computation is
needed. The same argument applies at the outer face by reflection.

For W=u(r)c with u=(r-a)^2 near the inner face and smoothly cut off before
the outer face, these defects are O(1) in the radial principal momentum term.
The smooth geometry, angular and lower-order terms contribute only O(h)
corrections to this defect because u=O(h^2) and d u/h=O(h) locally.
The normal STF principal symbol is c+er(er dot c)/3, which is invertible.
Consequently the failure is necessary for the actual continuum target too.
No additional traction or normal derivative is prescribed to mask it.

## Kernel and adjoint: logically separate conclusions

There is no moment-feasible resulting member whose kernel could be declared
absent. Kernel status for the requested consistent construction is UNRESOLVED
(nonexistent selection), not an assertion that every member has a kernel.
For illustration within this same family, setting the six free parameters
to zero detects the former vector (0,1,-1,0,...), yet the 16 by 14 scalar
Dirichlet derivative has exact rank 10. It therefore has four other scalar
kernel vectors. Each yields an exact full toroidal strain kernel W=u c
when angular-constant. This auxiliary rank check cannot rescue the failed
moments and is not a second candidate family or a positive feasibility claim.
No further search for parameters is performed.

For every formal member, form L from the original STF formula with these
radial derivatives. Let Hv=I3 tensor Hgeom, Ht=I9 tensor Hgeom,
Hgeom=diag(J)(h H tensor Htheta tensor Hphi). With positive angular weights,
these pairings are positive. Let Btraction=E^T Hb N pair vector traces with
original outward radial traction T n, n=-er,+er. Define
 Div T=Hv^-1[-L^T Ht T/2+Btraction T].
For arbitrary vector V and STF tensor T, multiplication gives exactly
 V^T Hv Div T=-(LV)^T Ht T/2+(EV)^T Hb(T n).
All nine tensor entries are counted. This is a full variable-coefficient
algebraic identity; it assumes no discrete geometry product rule.
For Dirichlet injection I, EI=0, M=I^T Hv I, L0=LI, hence
 M A=-L0^T Ht L0/2.
Adjoint identity PASS applies to the family even though consistency fails.
The original expanded production matrix is not asserted to satisfy it and is
not modified; the historical quadrature obstruction is respected.

Positive energy alone gives neither injectivity nor a mesh-independent Korn
bound. The illustrative rank-deficient member has no strict coercivity.
No consistent member exists in this family, so a plausible coercivity estimate
for the requested construction is UNRESOLVED. There is no positive result
addressing both prior obstructions. No physical solves or further family follow.
