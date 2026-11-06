# Run 43 — certified fixed-mesh mixed assembly

Base: `6b71b0f3913454f1cff92f01233b1d5102b1bccc`.
One agent. **Assembly CERTIFIED for the exact forms on this mesh.** No physical
loads or solves, nonlinear iterations, mesh sweeps, historical audit reruns,
Actions, or changes to the validated integration code or production solver.
This result makes no coupled-existence, positivity, or global stability claim.

## Smallest mesh and spaces

There is one parameter brick, the smallest possible nonempty partition:

    A(xi,eta,zeta)=(1/20+xi/5, 2*pi*eta, 2*pi*zeta).

Use the specification's exact X, with the two pairs of angular faces
self-identified and both radial faces physical. The aspect ratio is
10*pi < 32; simultaneous parameter bisection would retain this finite bound.
No refinement is executed. One periodic cell is allowed by the stated mesh
class; no extra aspect-ratio threshold or multiple-cells-per-period rule is
introduced. It is not a minimum-dimension theorem among every possible mesh.

The independent raw scalar generators are the 27 Q2 monomials in reference
coordinates and the nine nonconstant Cartesian monomials of degree at most
two. Cartesian constant 1 is already present. Orderings are in `result.json`.

Write b(t)=t(1-t). A basis of the entire periodic scalar space is

    xi^a b(eta)^j b(zeta)^k, a=0,1,2; j,k=0,1,
    x,y,z,x*x,x*y,x*z,y*y,y*z,z*z.

The first 12 functions impose the full polynomial periodic trace constraints;
the Cartesian enrichment is automatically periodic. Periodicity on a Q2
polynomial imposes equality of its linear and negative quadratic coefficients
in each angular variable, precisely giving span{1,b}. No nodal trace rule is
used. The four homogeneous functions are

    xi(1-xi) b(eta)^j b(zeta)^k, j,k=0,1.

Here is why these are the full constrained spaces, rather than a subspace
selected for convenience. A Cartesian quadratic pulled back to a radial
torus is a finite trigonometric polynomial. If it equals a polynomial in the
two angular coordinates on an open face, analytic continuation and
periodicity make that polynomial constant in each angle. At any positive
fixed radius, the nonconstant Cartesian quadratic coefficients are then all
zero. Explicitly, phi's second harmonics force equal xx,yy coefficients and
zero xy; its first harmonics remove xz,yz,x,y; theta's first and second
harmonics then remove xx=yy,zz,z. The executable independently checks rank
nine of the nonconstant Laurent/Fourier trace coefficient matrix at r=1/20.
Thus enrichment cannot cancel a nonconstant polynomial radial trace. It
vanishes from the homogeneous space; the Q2 radial factor must be xi(1-xi).
The same argument throughout the brick proves raw scalar independence.

`N` is the exact 36-by-21 periodic coefficient injection. `I` is the exact
21-by-4 homogeneous injection. `J` selects the remaining 17 basis functions
as boundary lifts. `[I J]` is invertible. Cartesian vector spaces use the
Kronecker product of these matrices with the 3-by-3 identity: scalar-major,
Cartesian-component-minor ordering. Dimensions are 63 and 12, respectively.

The tensor basis is the **same 98 functions**, with the same ordering, as
the preserved integration certificate, analytically extended over this
whole brick. In terms of its small-element reference coordinates y,

    y0=20*xi-5, y1=200*pi*eta, y2=200*pi*zeta.

An invertible affine coordinate substitution preserves Q2. Its differentiated
constant is zero, so the 78 nonconstant D generators span precisely the
specified Q2 strain-generator span. The remaining 20 E_beta*(1,x,y,z)
generators are unchanged. Independence follows already on the embedded
certified cell. There are no tensor face constraints (L2 conformity).

## Exact matrices, identities, and reactions

`L` (98-by-63) represents the exact full strain. Its Q2 part is constructed
by exact affine polynomial substitution and differentiation in the original
certified coordinates. The Cartesian quadratic part uses exact Cartesian
derivatives decomposed in the specified STF basis. All coefficient operations
are symbolic over rational numbers and pi. `L0=L*(I tensor I3)` has rank 12.
The program checks the Q2 substitutions, all 30 Cartesian vector polynomial
strains (including the three zero constant strains), entire traces, and the
radial trace rank. Polynomial reproduction concerns the full space, not
nonzero Cartesian quadratics satisfying homogeneous boundary data.

Let H denote the exact tensor mass, with full Frobenius contraction. The
authoritative assembled definitions are

    mixed = H*L,
    Dweak = -L0.T*H/2,
    C = L0.T*H*L0/2.

No displacement Riesz mass or strong divergence matrix is needed: Dweak is
the dual-space operator. These identities retain common exact H entries;
they are not independent interval entries that can be chosen arbitrarily.
The symmetric scaled mixed block is

    [ H          -H*L0 ]
    [ -L0.T*H       0  ].

Its Schur complement is -2C. This saddle matrix is symmetric and indefinite;
positive definiteness is certified for H and C, not for this saddle matrix
or for the full-space strain stiffness with its conformal Killing kernel.

The exact adjoint follows from T:grad(z)=T:Lz/2 for STF T. This is checked
symbolically for arbitrary entries. The executable checks the formal matrix
transpose identity and C symmetry. For any load functional f and tensor
coefficient sigma, the retained reaction vector is

    Jvec.T*(f+L.T*H*sigma/2).

Changing any lift by Ivec*c changes its value by
c.T*(Ivec.T*f+L0.T*H*sigma/2), exactly the interior equation residual. Thus
reaction lift independence is certified modulo the interior equation; no
physical right-hand side is chosen or solved. Reactions are not set to zero.
The reaction matrix is recoverable exactly from the stored mixed matrix and
J; the JSON records this definition. No numerical face-traction surrogate is
used for these specified variational reactions.

For completeness the executable also checks all three exact cofactor columns
of DX. This chart has det(DX)=-r*q. Physical outward oriented face measure is
sign(det B)*cof(B)*n_ref, while its magnitude is |cof(B)*n_ref|.
Consequently radial outward measures are -r*q*e_r and +r*q*e_r in dtheta dphi,
as required. Angular cofactors include the corresponding affine side factors.
Applying the divergence theorem to the pointwise identity yields the exact
adjoint including radial terms and paired seam jumps. The tensor basis need
not be periodic, so seam traction contributions must be paired as jumps,
not discarded. This is the specification's variational adjoint and reaction
check, not an independent numerical surface integration experiment.

## Integration uncertainty and positive definiteness

The preserved `run43_integration/result.json` is checked against its recorded
SHA-256. Its calculations are not rerun. For the new whole-shell H entries,
the executable reuses the certified tensor evaluator and the same validated
range-integration procedure: evaluate the entire mass integrand on one closed
reference box with 128-bit Arb arithmetic. Multiply its enclosure by the box
volume (one). Use |det B|=(4*pi^2/5)*r*q, with r>=1/20 and q>=3/4, and
Frobenius weights (1,1,1,2,2,2). All 9,604 mass entries, 6,174 mixed entries,
and 144 C entries are recorded with exact dyadic midpoint/radius bounds.
The mixed and C bounds follow by outward interval matrix arithmetic. Mirrored
C bounds use the upper entry, valid because the exact form is symmetric.
The common exact algebraic definitions remain authoritative.

These enclosures are deliberately coarse. They certify exact integral
containment, not a quadrature replacement, relative accuracy, or midpoint
positive definiteness. No midpoint eigenvalues are used. There is no omitted
integration remainder: the range bound includes every point of the box and
Arb accounts for arithmetic and transcendental rounding.

A rigorous structural certificate avoids unnecessary bound refinement. The
old certified physical cell K0 is contained in this single whole-shell cell,
and the tensor functions and coefficient ordering are identical. Therefore
for every real coefficient vector c,

    c.T H c = integral_Omega |T_c|^2
             >= integral_K0 |T_c|^2 >= 2^-332 ||c||^2.

This reuses the existing exact-mass certificate, including its integration
uncertainty, instead of asserting positivity for every arbitrary matrix in
the coarse entry boxes. The exact global mass is distinguished from those
arbitrary matrices. No small-cell eigenbound is asserted for a different
coefficient basis without transformation.

At 128 bits, a validated inverse of the exact coefficient Gram L0.T*L0
encloses its inverse. The sum of absolute inverse entries is strictly below
2^44, bounding its spectral norm. Thus

    lambda_min(L0.T*L0) > 2^-44,
    lambda_min(C) >= 2^-333 lambda_min(L0.T*L0) > 2^-377.

These are rigorous certificates for the exact forms on this one fixed mesh,
consistent with every certified integration uncertainty. Neither the stored
midpoint matrices nor arbitrary symmetric matrices within their entry boxes
are certified positive definite. No interval refinement was needed. No
mesh-uniform coefficient bound, nonlinear stability, or coupled solve follows.

## Reproduce and evidence

From repository root:

```sh
python -m pip install -r research/run43_assembly/requirements.txt
python research/run43_assembly/certify.py
```

`result.json` records the mesh, basis order, constraints and lift injections,
exact strain maps, exact operator definitions, every matrix enclosure and
error bound, inherited evidence digest, certificate exponents, and successful
assertions. Dependencies are python-flint 0.8.0 and SymPy 1.14.0. The tested
runtime is Python 3.12. Reproduction requires no network after dependencies
are installed. No preexisting file is modified.

All 18 executable checks passed. No required exact check failed. Independent
numerical surface quadrature and midpoint positive definiteness are not
claimed checks. The former is unnecessary for the specified variational
reaction definition; the latter is not required for the certified exact form.

One recommended next step: review this fixed-mesh certificate and its coarse
entry enclosures before authorizing any subsequent use. **STOP: await explicit
authorization before any coupled physical solve.**
