# Stage G: coupled stationary equations

2026-11-10 UTC. One Codex agent; bounded analytical stage. Branch
`research/stage-g-coupled-stationary-spectral`. Immutable parent:
`767783984a8502d892cf52048cd43eadde79e552`; baseline tree
`0da0d813d76108bc880d2e34c47f82b261026cd1` (1,120 tracked blobs).
The checkout was clean before branch creation directly at this parent. The
clone's default main tip was different; it was not used as the research parent.

**VERIFIED:** the exact reduction below. **CONDITIONAL / OBSTRUCTION IDENTIFIED:**
the companion report excludes nonzero real-frequency modes in a controlled
spherical-exterior, globally elliptic class. **NOT YET ESTABLISHED:** a general
coupled stationary exclusion or a localized solution. No theory changes.

## Frozen authority and notation

Preliminary inspection: original frozen
`evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`, equations
(3)–(9), (14)–(15); both Stage D reports, both Stage E reports, both Stage F
reports and the complete Stage F verifier. Stage D uses inherited static
fields; Stage E's additive phases do not exhaust Stage F rotations. The frozen
provenance verifier supports `--package-only`; its historical full-tree guard
is not applicable to later research additions. Recovery remains INCOMPLETE.
Every historical file is preserved; no earlier conclusion is silently extended.

With signature (-+++), G, lambda_V, v positive, the unchanged action is

    S = integral sqrt(-g)[R/(16 pi G) - (d rho)^2/2
                          - f(rho) d_mu Z* d^mu Z/2 - V(rho)],
    Z=(Theta-Theta_c)+i(Psi-Psi_c), f=1+a rho²>0,
    V=lambda_V(rho²-v²)²/4.

Complex notation describes precisely two real original fields; there is no
extra factor of two, mass, gauge field or added coupling. The phase chart must
support the rotation (real plane, or a suitable local image as in Stage F);
no global torus rotation is presumed. V'=lambda_V rho(rho²-v²), f'=2a rho.
Varying before reduction gives

    nabla_mu(f d^mu Z)=0,
    Box rho - V' - f' d_mu Z* d^mu Z/2=0,
    T_mu nu=C_mu nu-g_mu nu(C/2+V),
    C_mu nu=rho_mu rho_nu+f Re(d_mu Z* d_nu Z),
    C=g^mu nu C_mu nu,
    R_mu nu=8 pi G(C_mu nu+V g_mu nu).

## General stationary ADM reduction

No axisymmetry, circularity, staticity or restricted shift is imposed here.
All coefficients depend on x^i only, and

    ds²=-N²dt²+h_ij(dx^i+beta^i dt)(dx^j+beta^j dt),
    g^tt=-1/N², g^ti=beta^i/N²,
    g^ij=P^ij=h^ij-beta^i beta^j/N²,
    sqrt(-g)=N sqrt(h)=:w, N>0, h positive,
    n=(partial_t-beta^i partial_i)/N,
    K_ij=(D_i beta_j+D_j beta_i)/(2N), K=D_i beta^i/N.

D is the Levi-Civita connection of h; D_i on a scalar is partial_i.
Set rho=rho(x), Z=exp(i Omega t) z(x), Omega real constant, and define

    b=beta^i D_i z, p=(i Omega z-b)/N,
    q=-beta^i D_i rho/N,
    Pi_Z=exp(i Omega t) p,
    Xrho=|D rho|_h², Xz=h^ij Re(D_i z* D_j z),
    U=q²+f |p|², X=Xrho+f Xz,
    C=X-U, L_m=(U-X)/2-V.

In particular stationary rho does NOT imply q=0 when the shift is general.
The exact scalar equations are

    L_Omega z := (1/w) partial_i(w f P^ij partial_j z)
       +2 i Omega f beta^i partial_i z/N²
       +[i Omega/w partial_i(w f beta^i/N²)+Omega² f/N²]z=0,  (G1)

    (1/w) partial_i(w P^ij partial_j rho)
        -V' -a rho (Xz-|p|²)=0.                            (G2)

The equivalent unexpanded current form of (G1) is

    partial_i[w f (h^ij partial_j z+beta^i(i Omega z-b)/N²)]
                 -i Omega w f(i Omega z-b)/N²=0.           (G3)

For clarity all first coefficient derivatives in the spatial second-order
operator are

    f P^ij partial_i partial_j z
      +[partial_i(f P^ij)+f P^ij partial_i log w]partial_j z,
    partial_i(f P^ij)=f' rho_i P^ij+f partial_i h^ij
       -f (partial_i beta^i)beta^j/N²
       -f beta^i partial_i beta^j/N²
       +2f beta^i beta^j partial_i N/N³,
    w^-1 partial_i(w f beta^i/N²)
       = f partial_i beta^i/N²+f' rho_i beta^i/N²
         +f beta^i partial_i log w/N²-2f beta^i partial_i N/N³.

The same P and w appear in (G2), without a temporal rho derivative term.
The normal phase kinetic expansion is

    |p|²=[Omega²|z|²+|b|²-2 Omega Im(z* b)]/N².

Flat h=delta, N=1, beta=0, constant rho gives exactly
`L_Omega z=f (Delta+Omega²)z`; f>0 can be divided out.

## Complete stress, energy and momentum

All independent coordinate components, including every shift contribution
through L_m, are

    C_tt=f Omega² |z|²,
    C_ti=f Omega Im(z* D_i z),
    C_ij=rho_i rho_j+f Re(D_i z* D_j z),
    T_tt=f Omega²|z|²+(-N²+|beta|_h²)L_m,
    T_ti=f Omega Im(z* D_i z)+beta_i L_m,
    T_ij=C_ij+h_ij L_m.                                  (G4)

ADM quantities (j_i=-T(n,partial_i)) are

    epsilon=(U+X)/2+V,
    j_i=-q rho_i-f Re(p* D_i z)
       =[(beta^j rho_j)rho_i
          +f Re(b* D_i z)-f Omega Im(z* D_i z)]/N,
    S_ij=C_ij+h_ij[(U-X)/2-V],
    S=h^ij S_ij=(3U-X)/2-3V,
    epsilon+S=2U-2V,
    S_ij-h_ij(S-epsilon)/2=C_ij+V h_ij.                  (G5)

These are complete, without imposing j_i=0. Multiplication of all phase
values and derivatives by exp(i Omega t) leaves their Hermitian products
unchanged. Consequently every component in (G4) is stationary although
partial_t Z=i Omega Z need not vanish. This is a kinematic invariance proof;
no scalar or Einstein equation follows just from stationarity of T.

The exact positive slice matter energy is

    E_m=integral [q²/2+|D rho|²/2+V
                  + f(|i Omega z-beta.Dz|²/N²+|Dz|²)/2] dvol_h. (G6)

The last term is the rotational mode's positive energy. A vanishing Lorentz
contraction Xz-|p|² does not remove either positive contribution.

## Full stationary Einstein equations

With the above K convention the independent constraints are

    R[h]+K²-K_ij K^ij=16 pi G epsilon,                    (G7)
    D_j(K^j_i-delta^j_i K)=8 pi G j_i.                   (G8)

All six spatial evolution equations, set stationary, are

    0=-D_i D_j N+N[R_ij[h]+K K_ij-2K_ik K^k_j]
          +(Lie_beta K)_ij-8 pi G N(C_ij+V h_ij),        (G9)
    (Lie_beta K)_ij=beta^k D_k K_ij
                     +K_kj D_i beta^k+K_ik D_j beta^k.

The stationary trace evolution/lapse equation is

    D²N=N[K_ij K^ij+4 pi G(epsilon+S)]+beta^i D_i K
        =N K_ij K^ij+8 pi G N(U-V)+beta.DK.             (G10)

The lapse equation is redundant with the full system, but useful for integral
identities. In particular the beta.DK term is present when K is not constant.
Equivalently, Gauss–Codazzi and the normal Ricci identity give

    G_nn=(R[h]+K²-K_ij K^ij)/2,
    R_ni=-D_j(K^j_i-delta^j_i K),
    R_nn=D²N/N-K_ij K^ij-beta.DK/N=8 pi G(U-V).

The spatial spacetime Ricci projection is

    R4_ij=R_ij[h]+K K_ij-2K_ik K^k_j
                 -D_iD_j N/N+(Lie_beta K)_ij/N.

These formulas, the definition of K, both scalar equations, and (G7)–(G9)
constitute the coupled stationary system. The Hamiltonian constraint alone
is not equilibrium. Stage E is recovered by a purely axial circular shift
(K=0); Stage D is recovered only after beta=0 AND Omega z=0.

## Principal symbol and formal symmetry

For fixed real Omega, the principal symbol is f P^ij xi_i xi_j. In an h
orthonormal frame its eigenvalues divided by f are 1,1,
`1-|beta|_h²/N²`. It is positive elliptic iff the selected stationary Killing
field k=partial_t is timelike. It is degenerate at g_tt=0 and indefinite when
k is spacelike. N=0 or loss of h positivity would be a failure of the ADM
slicing; an ergosurface with N>0 is a different issue. The eigenvectors of
h and the coordinate location of an ergosurface can change under changes of
stationary generator. A rotating coordinate chart in Minkowski space can
make that generator spacelike: for that fixed generator the reduced principal
symbol really degenerates at its light cylinder, even though spacetime has
no physical gravitational ergoregion. Cross terms g_ti alone do not establish
an ergoregion. A stationary time redefinition t'=t+chi(x), with the same k,
changes z by exp(-i Omega chi); it does not remove g(k,k)=0 or its degeneracy.

On compactly supported complex profiles, L_Omega is formally Hermitian for
`<u,v>=integral w u* v d³x`. The divergence term is symmetric; with
`A^i=f beta^i/N²`, the shift term is `i Omega[2A.D+div_w A]`, also symmetric
because `2A.D+div_w A` is skew-adjoint. Omega² f/N² is real. This is formal
symmetry, not a positive quadratic form, a specified self-adjoint domain,
a standard frequency-independent eigenvalue problem, or a coercive elliptic
operator in an ergoregion. The frequency family is a quadratic operator pencil.
No reduction to a self-adjoint Helmholtz operator is assumed globally.

## Exact verification conventions

`python DEE_STAGE_G_VERIFIER.py` uses Python/SymPy exact arithmetic. Arbitrary
3-dimensional coefficient functions verify (G1) and (G2) by direct covariant
density differentiation; general pointwise frame jets verify the full stress
and projections. Block matrix factorization verifies general ADM inversion
and determinant. For curvature cost control, direct 4D Ricci is computed on
an exact stationary family with N,H,W,b,c arbitrary functions of x,

    h=diag(1,H²,W²), beta=(b,c,0).

This has a nonzero shift along coefficient gradients, generally K!=0, and
mixed K_xy; it is not the circular K=0 test of Stage E. Every Hamiltonian,
momentum, spatial Ricci and normal Ricci component is checked in that family.
It is a documented exact component reduction, not a claim that this family
represents every geometry. General equations also follow from Gauss–Codazzi
with stated signs. Symbolic checks certify local algebra only. ODE asymptotics,
unique continuation, positivity arguments and global hypotheses are analytical,
not certified by a PASS counter. Actual run records are in the companion file.
