# Stage I: exterior spectral theorem audit

2026-11-10 UTC. One agent; no delegation, parallel research, numerical PDE
solvers, simulations, or theory changes. Scope ends at Stage I.
Exact parent: `72f2ff834bde7d046c5104310aa8b0b9df0254ee`.
Branch: `research/stage-i-exterior-spectral-uniqueness`.

**Outcome C — exact unresolved proof/citation obligation.** This bounded audit
does not certify exterior uniqueness for the full H6 class. The coefficients
fit the asymptotic conditions of Koch–Tataru unusually closely. The remaining
obligation is a sign-consistent, parameter-uniform exterior Carleman estimate
for the full operator, specified below and in the companion report. This is
a limitation of this audit, **not a claim that the mathematical literature
has no applicable theorem, that H6 is insufficient, or that this spectral
problem is an open problem in the literature**. No discovery or counterexample
is claimed. The preprint's sign inconsistency does not invalidate its intended
positive-energy exclusion or supply Outcome D for Einstein–DEE.

## Immutable inputs and exact coefficient class

The frozen original source inspected was
`evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`, conventions
and equations (3)–(9), (14)–(15). The Stage H inputs inspected were both Stage H
reports, H1–H13 and the spectral/function-space discussion, plus the verification
requirements. Previous stages were not rederived.

The action in the authorized repeated-eigenvalue sector is unchanged:

    S = integral sqrt(-g)[R/(16 pi G) - (d rho)^2/2
                          - f d_mu Z* d^mu Z/2 - V],
    f=1+a rho^2>0, V=lambda_V(rho^2-v^2)^2/4, f_v=1+a v^2>0.

Signature (-+++), positive target kinetic metric, rho=rho(x),
Z=exp(i Omega t)z(x), real Omega !=0. Complex notation packages the original
two real phase fields; it adds no field or coupling. The Einstein equations,
including their source and sign conventions, remain unchanged.

With w=N sqrt(h), P=h^{-1}-beta beta/N^2, the exact operator is

    P_Omega z = div(A grad z)+i Omega(2B.grad+div B)z+Omega^2 C0 z,
    A=w f P/f_v, B=w f beta/(f_v N^2), C0=w f/(f_v N^2).

A is real symmetric; B and C0 are real. Work in Cartesian L^2(dx), not a
Dirichlet realization on the artificial interface. H6 means, exactly as in
Stage H, k>=3, delta>0, p0>=1, s>3/2, smooth angular coefficients, and

    N=1+n1(hat x)/R+O_k(R^{-1-delta}),
    h_ij=delta_ij+h1_ij(hat x)/R+O_k(R^{-1-delta}),
    beta^i=b1^i(hat x)/R^{p0}+O_k(R^{-p0-delta}),
    eta=rho-v=O_k(R^{-s}).

O_k(R^{-q}) controls each Cartesian derivative alpha through order k by
C_alpha R^{-q-|alpha|}. Smooth angular data include all their sphere derivatives.
In particular, for every |alpha|<=k, the exact coefficients obey

    |partial^alpha(A-I)| + |partial^alpha(C0-1)| <= C_alpha R^{-1-|alpha|},
    |partial^alpha B| <= C_alpha R^{-p0-|alpha|}.

These follow by differentiating products, inverse matrices and positive
square roots in the exact definitions, on a sufficiently far tail. They do
not discard shift squares, angular derivatives, or amplitude terms.
At leading order, with n=n1/R, H=h1/R, T=tr H and
F=(f-f_v)/f_v=O_k(R^{-s}), Stage H gives

    A-I=(n+T/2+F)I-H+O_k(R^{-min(2,1+delta)}),
    C0-1=F+T/2-n+O_k(R^{-min(2,1+delta)}),
    B=b1/R^{p0}+O_k(R^{-p0-min(1,delta)}).

On that tail A is uniformly positive definite, C0>0, and all densities and
Euclidean/metric volumes are comparable. H6 is an asymptotic hypothesis; it
does not assert ellipticity through an arbitrary compact interior or ergoregion.
The requested class is z in H^1(E) intersect L^2(E), satisfying the equation
weakly. Its interface trace is any admissible H^{1/2} trace; no zero trace
or Sommerfeld condition is assumed. Local elliptic regularity supplies
H^2_loc on the elliptic tail: bounded first-order coefficients make the
right side of the divergence equation L^2_loc, and Lipschitz uniformly
elliptic A permits the local difference-quotient estimate. No weighted
regularity or decay improvement follows merely from that statement.

## Actual literature read, and retrieval boundary

H. Koch and D. Tataru, *Carleman estimates and absence of embedded
eigenvalues*, Communications in Mathematical Physics **267** (2006),
419–449, DOI [10.1007/s00220-006-0060-y](https://doi.org/10.1007/s00220-006-0060-y).
The author publication list independently confirms this citation:
<https://math.berkeley.edu/~tataru/research.html>.

The **actual full preprint**, not only its abstract, was read at
<https://arxiv.org/html/math-ph/0508052v1> and checked against the downloaded
TeX source at <https://arxiv.org/src/math-ph/0508052> (the source is a gzip
stream containing em4.tex). The relevant inspected parts are Definitions 1–2,
Assumptions A1–A4, Theorems 3, 8, 11, Corollaries 9, 12, Proposition 4,
equations (3)–(14), (35)–(43), the proof of Theorem 3 Steps 1–3, and §7's
variable-coefficient proof discussion. These are preprint numbering; no
unseen published-page theorem numbering is asserted.

Retrieved source SHA256:
`60205e287dccb26d4791f5b7ca8ee405cd8cce09bf0e84996234a3b20d9929c6`.
Retrieved HTML SHA256:
`d24e3f87b41684176c036b95479845a2f3e8e63c94881f3d209c751b22a7d690`.
These identify the text inspected, not a frozen-package addition or a proof
of the paper's mathematics. Source material was kept outside the repository.

The unversioned arXiv PDF endpoint returned HTML; an attempted v2 endpoint
also returned HTML (the inspected submission has v1). The Springer PDF URL
returned HTML, not a verified published PDF. An initially guessed author PDF
and papers index returned 404; the actual research index was then used.
A bounded supplementary attempt to retrieve Aronszajn's classical local
continuation paper from Numdam timed out on a PDF and an article-form URL
returned HTML. No theorem statement from that failed retrieval is imported.
No subsequent literature search was conducted. These retrieval failures are
not regression-test results.

## Verified statements and mapping

Koch–Tataru use W^{s,p} with Bessel-potential norm and dyadic
l^2 W^{s,p} norms over overlapping annuli {2^{j-1}<=R<=2^{j+1}}.
For n=3 the exponents in (7) are W^{1/4,4} and W^{-1/4,4/3}.
Their C^2_{<x>} norm bounds |g|, <x>|Dg|, <x>^2|D^2g|.

**Theorem 3:** for (-Delta-V)u=Wu on all of R^n, A1 requires V in
C^2_{<x>}, liminf V>0 and

    tau0=-liminf_{R->infinity} (x.grad V)/(4V) < 1/2.

A2 splits W=W1+W2, with small limiting dyadic multiplier norm X for W1 and
small limsup R|W2| for W2. X is the multiplier norm
W^{1/(n+1),2(n+1)/(n-1)} to W^{-1/(n+1),2(n+1)/(n+3)}.
For sufficiently small bounds and tau1>tau0, the stated weighted H^1_loc
solution is zero. Theorem 3's displayed weight is
(1+R^2)^{tau1-1/2}u in L^2; Proposition 4 instead uses
R^{tau1-1/2}v in L^2. We do not silently identify those two powers.
Choosing tau1=1/4 below makes unweighted L^2 sufficient for both.
This scalar-principal theorem alone does not handle div((A-I)grad).
The basic L^2-potential result also cannot absorb a nonzero Coulomb term
in n=3: R^{-1} is not in L^2(E). That fact does **not** rule out the paper's
long-range or variable-coefficient extensions.

**Theorems 8 and 11:** the additional gradient potentials Z^l grad u and
div(Z^r u) satisfy small limiting dyadic L^n bounds. A3 requires a low-frequency
part S_{<M}Z to satisfy A2; A4 replaces that requirement by
<D>^{-N}Z satisfying A2 for sufficiently large N. They are not assumptions
on an arbitrary imaginary scalar potential. For the mapping below these
extra potentials are zero, so their multiplier, integrability and smoothing
conditions are satisfied identically. In particular the p0=1 shift is not
incorrectly forced into the short-range Z class.

**Section 7, equation (40):** a real-coefficient elliptic self-adjoint
operator is written

    L=-partial_j(a^{jk}partial_k)+i(b^j partial_j+partial_j b^j)+c.

The first-order expression is the Hermitian sum of a derivative and its
product-composition partner: i(2b.grad+div b). The assumptions are

    a,b,c in C^2_{<x>},
    limsup R|grad a| <= delta0,
    limsup (|b|+R|grad b|) <= delta1,
    limsup (|c|+R|grad c|) <= delta1^2,

with sufficiently small constants. Equation (43) additionally normalizes
a to be close to I at infinity; our a already tends to I and needs no
coordinate transformation.

For fixed real Omega, the **exact candidate dictionary** is

| Literature quantity | Stage H quantity | Checked consequence of H6 |
| --- | --- | --- |
| a^{ij} | A^{ij} | real symmetric, elliptic on a tail, C^2_{<x>}, A->I, R grad A->0 |
| b^i | -Omega B^i | real, C^2_{<x>}, |b|+R|grad b|->0 even for p0=1 |
| c | Omega^2(1-C0) | real, C^2_{<x>}, |c|+R|grad c|->0, permits Coulomb 1/R |
| L | H_Omega | exact Stage H fixed-frequency operator |
| positive spectral parameter | Omega^2 | fixed and strictly positive |
| long-range V in (-Delta-V) convention | Omega^2 | A1, tau0=0 |
| W, Z | 0, 0 | A2 and A4 identically satisfied |
| tau1 | 1/4 | tau1>0; both displayed weighted L^2 requirements follow from L^2 |

Here

    H_Omega=-div(A grad)-i Omega(2B.grad+div B)+Omega^2(1-C0),
    Q_Omega=H_Omega-Omega^2=-P_Omega.

Consequently all limiting constants in (40) can be zero; equivalently any
specified positive small constants hold eventually. Large fixed Omega only
changes how far out this occurs. H6 supplies the second derivative bounds,
not just pointwise convergence. Arbitrary angular tensor data are retained.
Coulomb terms are placed in c, which (40) permits; alternatively their smooth
long-range scalar nature is compatible with A1. They are not claimed small
in an L^2 potential norm. Formal Hermiticity supplies the displayed algebra,
not spectral positivity or a physical boundary domain.

**Corollary 12:** for a global self-adjoint elliptic realization satisfying
(40), sufficiently small delta0 and W,Z as in A2,A4 with delta=0, the text
states that P+W has no eigenvalues lambda>C delta1. Thus our dictionary is
strong positive evidence for absence of global eigenvalues of a suitable
globally elliptic extension. It is not itself a theorem about an exterior
solution with arbitrary interface trace or an inhomogeneous cutoff equation.
No assumption about such a global physical realization is established here.

## Exact source-sign issue and the remaining applicability obligation

The v1 TeX as well as HTML define L with **negative** divergence principal
part but print equation (42), in Theorem 11, as

    Lu+Vu=Wu+Z^l grad u+div(Z^r u).

The physical dictionary instead gives (L-Omega^2)z=0. Taking V=-Omega^2
violates A1's positive lower limit. Taking V=+Omega^2 in the printed (42)
gives a different equation. Taking c=-Omega^2 C0 would give a nonzero
limiting c and lose the vanishing-small-constant mapping just verified.
The earlier equation (3) is (-Delta-V)u=Wu, and Corollary 12 concerns
positive eigenvalues, indicating an intended opposite sign in (42).
**ERROR IDENTIFIED:** this internal sign inconsistency in the inspected
preprint; no assertion is made about the uninspected published equation.
We do not silently edit a cited theorem and call its application verified.

The proof in §7 says, verbatim, “Once (43) is gained the Carleman estimates
(7) in Proposition 4 remain valid with essentially no change.” It also says
compact-region geometry may require extra convexity of the weight. This is
important positive evidence, and may be enough for a specialist's completed
application. This audit has not verified the sign-correct operator version
and its parameter-uniform exterior support/weight formulation through those
modifications. In particular, it does not promote this prose to a newly
proved uniform estimate or treat Corollary 12 as that estimate. The unresolved
obligation is narrower than a new general scattering theory: certify the
specific inequality (I*) in the companion report for Q_Omega, with the exact
solution class, no boundary condition at R0, and constants independent of
the large Carleman parameter and taper parameter. A sign-consistent checked
published statement or a rigorous verification of this extension would
potentially settle this obligation; neither is claimed accomplished here.

Flat exterior Rellich is limited to the exact Helmholtz equation and is not
an applicable replacement for this tensor operator. The classical scalar
criteria and §7 were examined within this single source rather than opening
an unrelated search. No flat, radial, outgoing-only, Dirichlet, positive-
commutator-modulo-compact-errors, or abstract self-adjointness assertion is
used to claim the requested exclusion.

## Phase 1 stopping decision

**VERIFIED:** full source read, exact coefficient mapping and differentiated
bounds above. **CONDITIONAL:** the intended variable-coefficient positive-
energy estimate would fit H6 on a tail. **NOT ESTABLISHED:** a fully certified
exterior theorem with arbitrary trace. Proceeded once to the cutoff Carleman
attempt, then stopped at the precise estimate obligation. The conclusion is
an audit limitation, not a negative result about Koch–Tataru or a proof that
anisotropy creates embedded modes.

The Einstein equations have not been shown to imply H6 from finite ADM mass
and finite matter energy. Stage H's leading Einstein constraints are
conditional on H6 and do not supply its full differentiated bootstrap. Even
successful spectral exclusion within H6 would therefore remain conditional
as a claim about all finite-energy Einstein–DEE configurations. Continuation
through an interior ergoregion is a separate issue and is not asserted.
