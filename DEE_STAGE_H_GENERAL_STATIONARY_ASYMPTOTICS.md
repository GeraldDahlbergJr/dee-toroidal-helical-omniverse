# Stage H: general stationary asymptotics

2026-11-10 UTC. One Codex agent, no subagents or parallel investigations.
Exact parent: `ca30828547b5989c751b5c9639796377cdc1e3c9` (Stage G).
Research branch: `research/stage-h-general-stationary-spectral`.
Baseline tree: `3c7f9c73bb4a3134ba84840452cde5020c4c0e28`;
**1,123 tracked paths**, directly counted, not inferred from earlier reports.
The complete NUL-delimited recursive tree (path, mode, type, object ID) was
recorded outside the repository and is reproducible with
`git ls-tree -rz ca30828547b5989c751b5c9639796377cdc1e3c9`.
The immutable Git tree is the complete baseline inventory, including all
archives; a selective inventory is not substituted for it.

**Scientific classification: Outcome C — necessary conditions and residual
spectral class.** **VERIFIED:** the general density pencil, first variations,
conditional leading Einstein constraints and noncircular Komar cancellation
below. **CONDITIONAL:** differentiated end expansions and mass equality.
**NOT YET ESTABLISHED:** a general anisotropic positive-frequency exclusion,
or a localized solution. The companion spectral report explains why formal
symmetry and the Coulomb logarithm alone do not prove that exclusion.

## Preliminary authority and immutable theory

The origin URL was verified as
`https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse.git`.
The Stage G commit exists as a commit object. A new local branch was created
with `git switch -c research/stage-h-general-stationary-spectral
ca30828547b5989c751b5c9639796377cdc1e3c9` (one command, shown wrapped).
The default main tip was different and was not used as the parent. The
requested remote Stage H branch did not exist at preliminary inspection.
No AGENTS.md was found in the checkout or its workspace parent.

Before new derivations, inspection covered the frozen original
`evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`, especially
(3)–(9), (14)–(15); Stage D's two reports (static phase and gravitational
obstructions); Stage E's two reports (additive stationary representation,
circularity and rotation boundary hypotheses); Stage F's two reports
(target Killing algebra, fixed vacuum and rotational residual); both Stage G
reports and the complete Stage G verifier. Equations G1–G10 are inputs checked
again below where used; their hypotheses are not silently strengthened.
The frozen verifier was inspected: its supported later-stage mode is
`python evidence/einstein_dee_sources/verify_evidence.py --package-only`.
It checks the frozen package but explicitly skips Git preservation. The
separate exhaustive baseline comparison in the companion report supplies that
check. Recovery remains INCOMPLETE; physical tests remain NOT RUN.

Signature (-+++), G>0, lambda_V>0, v>0. In the authorized A=aI sector the
unchanged action and equations are

    S = integral sqrt(-g)[R/(16 pi G) - (d rho)^2/2
                           - f d_mu Z* d^mu Z/2 - V],
    f=1+a rho²>0, V=lambda_V(rho²-v²)²/4,
    Z=(Theta-Theta_c)+i(Psi-Psi_c),
    nabla_mu(f d^mu Z)=0,
    Box rho - lambda_V rho(rho²-v²) - a rho d_mu Z* d^mu Z=0,
    C_mu nu=rho_mu rho_nu+f Re(d_mu Z* d_nu Z),
    T_mu nu=C_mu nu-g_mu nu(C/2+V),
    R_mu nu=8 pi G(C_mu nu+V g_mu nu).                    (H1)

Complex notation is exactly the original two real phases. There is no
additional field, mass, potential or kinetic modification. Rotations require
a real phase chart or admissible local image as in Stage F; no global rotation
of an arbitrary phase torus is asserted. V gives amplitude mass squared
mu²=2 lambda_V v², while the phase Hessian is zero at the vacuum.

## Exact ADM reduction and density pencil

No axisymmetry, circularity, staticity, conformal flatness, prescribed
multipoles, or global timelike Killing field is imposed. Set

    ds²=-N²dt²+h_ij(dx^i+beta^i dt)(dx^j+beta^j dt),
    g^tt=-N^-2, g^ti=beta^i/N²,
    g^ij=P^ij=h^ij-beta^i beta^j/N², w=N sqrt(h),
    rho=rho(x), Z=exp(i Omega t)z(x), Omega real !=0.

Direct differentiation of `partial_mu(w f g^mu nu partial_nu Z)=0`
gives exactly

    L_Omega z = w^-1 partial_i(w f P^ij partial_j z)
       +2i Omega f beta^i partial_i z/N²
       +[i Omega w^-1 partial_i(w f beta^i/N²)
                                  +Omega² f/N²]z=0.     (H2)

This checks G1 independently from the original phase equation; in particular
both time-space cross derivatives contribute. Define f_v=1+a v²>0 and real
densities (these are not new physical fields)

    A^ij=w f P^ij/f_v, B^i=w f beta^i/(f_v N²),
    C0=w f/(f_v N²),
    mathcal P_Omega := (w/f_v)L_Omega
       =partial_i(A^ij partial_j)
            +i Omega(2B^i partial_i+partial_i B^i)+Omega² C0. (H3)

C0 here is a scalar density coefficient, distinct from H1's kinetic trace C.
The exact flat limit is `(Delta+Omega²)z=0`. The tensor principal quadratic
form of L is `f h^ij xi_i xi_j-f(beta.xi)²/N²`. In an h-orthonormal frame
its eigenvalues/f are `1,1,1-|beta|_h²/N²`. Thus it is elliptic precisely
where k=partial_t is timelike, degenerates at its ergosurface, and has
indefinite principal form where k is spacelike. This statement is about the
fixed normalized generator, not merely whether a coordinate cross term exists.
N>0 and h>0 need not fail at an ergosurface.

For arbitrary complex u,v and real Omega the local Green identity is

    w[u* L v-(L u)* v]
      =partial_i{w f P^ij(u* partial_j v-v partial_j u*)
                      +2i Omega w f beta^i u* v/N²}.    (H4)

It proves formal Hermiticity on compact supports in L²(w dx). H3 is instead
formally Hermitian in L²(dx): multiplication by w/f_v exactly supplies the
density. Dividing H3 by another spatial function need not preserve that
measure. Formal Hermiticity does not choose a self-adjoint domain and does
not imply positivity; the family is a quadratic pencil.

## What the basic assumptions actually supply

Assume one normalized asymptotically Euclidean end, N->1, h->delta,
beta->0 in an inertial end chart, rho->v and z->0, smooth stationary coupled
equations, finite ADM mass and finite proper matter energy. The limits and
ADM mass require a specified admissible asymptotic chart and first-derivative
surface convergence; they do not specify all higher derivative bounds.
For geometric arguments requiring completeness, no horizons or physical
inner boundaries, those are additional global hypotheses, not consequences
of an end alone. DEC holds for the positive target kinetic metric and V>=0.
For a unit observer the energy is the sum of half the squared normal and
spatial target derivatives plus V; the associated flux obeys DEC by
Cauchy–Schwarz (and sums of canonical real-scalar stresses preserve DEC).

Let

    p=(i Omega z-beta.Dz)/N, q=-beta.Drho/N,
    U=q²+f|p|², X=|Drho|²+f|Dz|²,
    epsilon=(U+X)/2+V, j_i=-q rho_i-f Re(p* D_i z).

Finite E_m=integral epsilon dvol_h implies Dz in L², p in L² on the end.
The identity `i Omega z=Np+beta.Dz` gives

    Omega² |z|² <= (N²+|beta|_h²)(|p|²+|Dz|²).

Bounded N,beta and f bounded below give z in H¹ intersect L². Near the
vacuum V is comparable to (rho-v)², so eta=rho-v and Deta are in L² there.
These are **derived integrated estimates**, not pointwise differentiated
power laws. DEC bounds orthonormal stress by a multiple of epsilon, so stress
is integrable in that frame; it does not turn integrable spikes into uniform
O_k(R^-3-delta) decay. Smoothness at every point supplies no uniform end
bounds on all jets. Finite ADM mass alone supplies neither a complete 1/R
series nor an angular-independent leading coefficient in arbitrary charts.

The exact forced amplitude equation is

    [w^-1 partial_i(w P^ij partial_j)-mu²]eta
      =lambda_V(3v eta²+eta³)+a(v+eta)(|Dz|²-|p|²).       (H5)

Only L¹ control of this quadratic source follows directly from the energy
estimate with bounded end coefficients. Claiming a Yukawa tail, uniform
bounds on derivatives of f, or a full stationary metric expansion from it
would require a further weighted elliptic regularity/bootstrap argument.
No such unproved bootstrap is used here. This identifies a proof gap, not
a counterexample to the coupled equations enforcing better decay.

## A controlled anisotropic class and its constraints

Write O_k(R^-s) for Cartesian derivative bounds through k with each
spatial derivative gaining one power. For expansion calculations only,
assume k>=3, delta>0, p0>=1 and s>3/2, with smooth angular functions,

    n=n1(hat x)/R, H_ij=h1_ij(hat x)/R,
    b^i=b1^i(hat x)/R^p0,
    N=1+n+O_k(R^-1-delta),
    h_ij=delta_ij+H_ij+O_k(R^-1-delta),
    beta^i=b^i+O_k(R^-p0-delta), eta=O_k(R^-s).         (H6)

These are **assumptions**, including the exponents, rather than a theorem
from finite ADM mass. The letter H_ij in H6 is a perturbation, not the radial
conformal factor in the Stage G specialization. For p0<1 the exact beta beta
term in P can compete with or dominate 1/R; no version of the expansion below
silently drops it. H2–H3 remain valid without H6.

With T=delta^ij H_ij and d=n+T/2, let
F=(f-f_v)/f_v=(2av eta+a eta²)/f_v. The first density jets are

    sqrt(h)=1+T/2+O_k(R^-q0), w=1+d+O_k(R^-q0),
    P^ij=delta^ij-H^ij+O_k(R^-q0)+O_k(R^-2p0),
    A^ij=delta^ij+(d+F)delta^ij-H^ij+E_A^ij,
    C0=1+F+T/2-n+E_C,
    B^i=b^i+E_B^i,                                      (H7)
    q0=min(2,1+delta),
    E_A,E_C=O_k(R^-q0), E_B=O_k(R^-p0-min(1,delta)).

The E bounds include products with F and shift squares since s>3/2,p0>=1;
F itself is kept exactly. Tensor derivatives of the remainders in H3 gain
one power, so their first-order divergence coefficients are respectively
O_{k-1}(R^-q0-1) and O_{k-1}(R^-p0-min(1,delta)-1).
Substitution into H3 gives the anisotropic equation

    0=(Delta+Omega²)z
       +partial_i([(d delta^ij-H^ij)] partial_j z)
       +Omega²(T/2-n)z
       +partial_i(F partial^i z)+Omega² F z
       +i Omega[2b^i partial_i+partial_i b^i]z
       +partial_i(E_A^ij partial_j z)
       +Omega² E_C z+i Omega[2E_B^i partial_i+partial_i E_B^i]z. (H8)

This is a long-range *differential* perturbation. Its O(1/R) tensor times
second derivatives cannot be replaced with a scalar potential. If one instead
divides H2 by f, its exact amplitude dependence is `grad log f . P grad`
and its associated shift divergence term; the F Helmholtz-density term in
H8 partly represents multiplication by f. At first order the metric part of
H2/f is

    -H^ij partial_i partial_j z
      +[partial^j(n+T/2)-partial_i H^ij]partial_j z
      -2 Omega² n z,

plus `grad F . grad z` and shift terms. H8 and this expression differ by
`(d+F)(Delta+Omega²)z` at first order, as checked symbolically. An apparent
amplitude-dependent phase mass in one density normalization is not a new
vacuum mass.

The unchanged Einstein equations restrict the permitted angular data. With
K_ij=(D_i beta_j+D_j beta_i)/(2N), K=D_i beta^i/N, they are

    R[h]+K²-|K|²=16 pi G epsilon,
    D_j(K^j_i-delta^j_i K)=8 pi G j_i,
    0=-D_iD_j N+N(R_ij[h]+K K_ij-2K_ik K^k_j)
                         +(Lie_beta K)_ij-8 pi G N(C_ij+Vh_ij),
    D²N=N|K|²+8 pi G N(U-V)+beta.DK.                    (H9)

The last line is the trace evolution equation; deleting beta.DK in a general
slicing is incorrect. H9 has the G7–G10 signs with j=-T(n,.). It follows from
the same H1 trace reversal and Gauss–Codazzi conventions; the unchanged Stage G
verifier checks their curvature signs on its stated noncircular exact family.
We use no extrapolation of that family to prove end regularity.

Here is a precise leading-order consequence **conditional on H6**. H9 and
p0>=1 imply every orthonormal stress component has an expansion at order
R^-3 with uniform o(R^-3) remainder: G_nn and G_ni use two derivatives of
h,beta, and the spatial equations and lapse use two derivatives of N,h.
If p0>1 the momentum leading term is smaller than R^-3; at p0=1 it is of that
order. The coefficient of epsilon at R^-3 must be nonnegative. Finite matter
energy forces that coefficient to vanish: any positive smooth angular patch
would give a divergent logarithmic radial integral. DEC then forces all stress
coefficients at that order to vanish, including the trace-reversed source.
Thus within this controlled class the needed stress=o(R^-3) follows from
H9, integrability and DEC; it is not imposed independently there. This
argument does not supply H6 itself, nor improved matter derivatives.

Writing Rlin_ij(H) for

    [partial_k partial_i H_kj+partial_k partial_j H_ki
                         -Delta H_ij-partial_i partial_j T]/2,

the order R^-3 equations are

    Delta n=0,
    Rlin_ij(H)=partial_i partial_j n,
    partial_i partial_j H_ij-Delta T=0.                  (H10)

K², beta.DK and Lie_beta K are O(R^-4) or faster for p0>=1.
Since `Delta[n1(hat x)/R]=R^-3 Delta_S² n1`, smoothness on the entire sphere
gives n1 constant. Write it as -m_N. The h1 functions cannot be arbitrary:
they must solve the tensor equation in H10, in addition to the ADM surface
conditions. A solution is H_ij=2m_N delta_ij/R. General homogeneous solutions
differ by a linearized Ricci-flat tensor: in three dimensions linearized
Ricci-flatness implies vanishing linearized Riemann. Locally this is a pure
coordinate strain `partial_i xi_j+partial_j xi_i`; a global admissible gauge
with the required bounds and any parity/charge normalization still needs to
be chosen and checked. We do not assume such gauges preserve every allowed
ADM definition. This local gauge statement does not assert an all-orders
stationary expansion or a global chart construction.

For p0=1 the momentum leading equation is

    Delta b_i-partial_i partial_j b^j=0                 (H11)

at order R^-3. For p0>1, finite energy only gives j=o(R^-3), which is too
weak to discard a source of order R^-p0-2. In that case the leading equation is
`(Delta b_i-grad_i div b)/2=8 pi G j_i` at that order **only if j admits the
corresponding expansion**; a homogeneous H11 requires the extra assumption
j=o(R^-p0-2). In particular the spin tail is not derived from finite energy.
The curl-free time-slicing part b=grad chi solves the homogeneous equation.
For explicit angular constraints define

    D_i^(q) F = (delta_ij-hat x_i hat x_j)partial_hatx_j F
                                      -q hat x_i F,
    partial_i(R^-q F)=R^-q-1 D_i^(q)F.

Then H10 is a linear angular tensor system using
D_i^(2)D_j^(1), and H11 uses D_i^(p0+1)D_j^(p0).
These formulas specify the constraints without declaring arbitrary anisotropic
coefficients or decay exponents to be admissible Einstein data.

## Gauge, rest frame, mass monopole and spin

A stationary spatial change x'=x+xi changes the first metric jet by
`H'_ij=H_ij-partial_i xi_j-partial_j xi_i` in the same asymptotic frame.
Bounded angular xi can change 1/R anisotropy. Translations affect subleading
terms; rotations rotate the tensor and angular coefficients. Logarithmic
coordinate changes and angle-dependent time slices need separate admissibility
checks. Parity conditions are especially relevant for center-of-mass and
angular momentum; finite ADM energy does not supply those charges automatically.

For t'=t+chi(x), k is unchanged and

    z'=exp(-i Omega chi)z,
    g'^t'i=g^ti+P^ij chi_j,
    g'^t't'=g^tt+2g^ti chi_i+P^ij chi_i chi_j.            (H12)

The density wave equation is covariant, so these transformations transport
H2 exactly on spacelike new slices. To first order about the inertial metric,
beta'=beta+grad chi; its induced profile phase is essential. This does not
remove the sign of g(k,k), nor a non-exact gravitomagnetic one-form. Setting
beta=0 is not a free gauge choice for a general stationary metric. A boost
mixes t with x and generally destroys stationarity relative to the original
coordinate generator for a localized configuration; it is not an unrestricted
stationary gauge transformation. Normalized stationary rest-frame end,
vanishing ADM linear momentum, admissible gauge and sufficient decay must be
proved or assumed rather than inferred from the word stationary alone.

If an admissible rest gauge additionally has

    N=1-m/R+O_k(R^-1-delta),
    h_ij=(1+2m/R)delta_ij+O_k(R^-1-delta),
    beta=O_k(R^-2),

the ADM surface formula gives M_ADM=m/G; lapse Komar flux gives m/G as well.
H7 becomes

    A=I+F I+O_k(R^-q0),
    C0=1+4m/R+F+O_k(R^-q0), B=beta+O_k(R^-3),           (H13)

with q0=min(2,1+delta). The equation has `+4m Omega²/R`; the corresponding
Schrodinger effective potential has `-4m Omega²/R`. This coefficient is
**CONDITIONAL on normalization**, not a universal Cartesian tensor formula.
For example the areal Cartesian Schwarzschild jet is
H_ij=2m hat x_i hat x_j/R, n=-m/R: its density frequency coefficient is
2m/R and its radial principal correction is -2m/R. Dividing the radial
principal coefficient recovers optical `1+4m/R`; ignoring it would give a
wrong factor of two. The two spatial jets differ by the strain of xi=m hat x,
verified exactly. Schwarzschild is a coefficient/gauge comparison, not a
regular DEE soliton. The invariant content is long-range gravitational optical
propagation (and its logarithmic time delay) relative to the normalized k;
a density's scalar 1/R term alone is not invariant.

With an independently justified stationary spin gauge,
`beta=-2G(J cross x)/R³+O_k(R^-3)`, the leading shift operator is

    2i Omega beta.grad=-4i G Omega (J cross x).grad/R³.

For standard quantum angular momentum L=-i x cross grad this equals
`+4G Omega J.L/R³`. For J along the z axis and exp(i l phi) it is
`+4G Omega J l/R³`. This fixes the convention explicitly; no axial symmetry
is used in H2 or H8. A general b1/R^p0 instead gives a first-order
O(R^-p0) coefficient and divergence O(R^-p0-1), with no universal spin
interpretation. Shift-squared principal terms start at O(R^-2p0) in H6.
The mass monopole, tensor perturbations, spin derivative and amplitude
derivative corrections play distinct spectral roles.

## ADM/Komar reconciliation without circularity

This derivation starts with the exact Komar current, not a circular metric.
For the Killing field k normalized to unit time at infinity set
F^mu nu=nabla^mu k^nu (antisymmetric). The Killing identity is
`nabla_mu F^mu nu=-R^nu_lambda k^lambda`. With the standard oriented
binormal and sign calibrated to a positive Schwarzschild mass, define

    M_K[k;S]=-(8 pi G)^-1 integral_S nabla^mu k^nu dS_mu nu
            =(4 pi G)^-1 integral_S (D_s N-K_sj beta^j)dA. (H14)

Here s points toward the asymptotic end on an outer sphere. For an inner
surface define M_K with s pointing from the excluded interior toward the
exterior, so the domain boundary orientation is its negative. The ADM
formula uses our K=(-partial_t h+Lie_beta h)/(2N) convention. The normal
acceleration gives D_s N; the normal/tangential derivative of k supplies
`-K_sj beta^j`. Neither maximality nor circularity is used. Integrating the
Killing identity, using Einstein H1, gives on every compact domain D

    M_K[k;S_out]-sum M_K[k;S_inner]
      =2 integral_D (T_mu nu-T g_mu nu/2)n^mu k^nu dvol_h
      =integral_D [N(epsilon+S)-2 beta^i j_i]dvol_h.     (H15)

The factor two and the sign of j follow from n.k=-N and
T(n,beta)=-beta.j. Equivalently the trace-reversed source in H1 gives
`2[N(U-V)-beta.j]`. With rho stationary and the rotational representation,

    N U-beta.j
      =N q²+q beta.Drho
          +f[N|p|²+Re(p* beta.Dz)]
      =f Re[p*(Np+beta.Dz)]
      =Omega f Im(z* p).                               (H16)

The amplitude terms cancel because Nq=-beta.Drho; all shift contributions
are retained before cancellation. There is no assumption that p=i Omega z/N
or that beta is axial. The exact original Noether current and its charge
convention are

    J_rot^mu=f Im(Z* d^mu Z),
    n_mu J_rot^mu=f Im(z* p),
    Q_rot=integral n_mu J_rot^mu dvol_h,
    sqrt(-g) J_rot^t=-sqrt(h) f Im(z* p).              (H17)

This sign makes Q positive for beta=0,Omega>0. Other conventions reverse
both the charge definition and its displayed coefficient. Local conservation
holds on shell; conserved global charge also needs no flux through boundaries.
On a bounded normalized end H¹ ensures charge integrability. Taking the
exhaustion limit when all integrals and surface limits exist yields

    M_K[infinity]=sum M_K[inner]
                     +2 Omega Q_rot-2 integral N V dvol_h. (H18)

If other ends exist they carry their oriented Komar fluxes too. A regular
filled axis is not an inner boundary; artificial interfaces cancel. A regular
center has no limiting flux. A horizon is an actual omitted term. When there
is a Killing horizon with generator chi_H=k+Omega_H phi and the usual regular
axial normalization, linearity gives
`M_K[k;H]=kappa A_H/(4 pi G)+2 Omega_H J_H` with
`J_H=(16 pi G)^-1 integral_H nabla^mu phi^nu dS_mu nu` in the matching
orientation. This is an illustration only where those horizon/axis hypotheses
hold; no such structures are assumed for general noncircular interiors.
At an outer boundary for a generator changed to k+omega phi, its Komar mass
changes by -2 omega J. For the fixed inertial k in H14 there is no extra
bulk spacetime angular momentum term to append to H16. Q_rot is internal
rotation charge and is not ADM angular momentum.

For a second check using H9, on a finite domain integration of the lapse gives

    integral_boundary D_s N
       =integral_D N(|K|²-K²+8 pi G(U-V))
                         +integral_boundary K beta.s,
    integral_D N(|K|²-K²)
       =integral_boundary beta^i(K^j_i-delta^j_i K)s_j
                          -8 pi G integral_D beta.j.

The boundary trace terms cancel, leaving precisely H14–H16; discarding
beta.DK would obstruct this cancellation. In a rest end beta=O_1(R^-p0),
p0>=1, K=O(R^-p0-1), so the shift Komar surface term is
O(R^(1-2p0))->0. A lapse flux alone is not the exact general finite-surface
Komar expression.

**Resolution:** the identity in the question,

    M_ADM=2 Omega Q_rot-2 integral N V dvol_h,

is **CONDITIONAL but genuinely noncircular**. It requires no inner Komar
flux (regular complete horizonless slice with no physical boundary suffices),
convergent integrals, and equality M_ADM=M_K[infinity] for the normalized
stationary end. An explicit sufficient end normalization is H13, together
with its differentiated decay and the ADM/Komar surface definitions: direct
surface evaluation then gives both m/G. More generally a stationary
ADM–Komar equality theorem requires its own admissible coordinates, decay,
field equations and boundary hypotheses. We do not invoke such a theorem
from finite ADM mass alone. Without those conditions the correct statement
is H18, plus the possible difference M_ADM-M_K[infinity].

Thus the narrower suggestion that noncircularity itself invalidates the
cancellation is refuted by H16. Stage G's G22 is correct with its expressly
stated boundary and end mass-equality qualifications. Its unrestricted
extension without those qualifications would be incorrect. No historical
file is edited. Matter-energy positivity and positive ADM mass do not give
spectral positivity: even beta=0 leaves the positive term
`2 Omega² integral f|z|²/N` balancing the nonnegative potential term.
There is no static-sign contradiction and no use of the positive-mass theorem
as a real-frequency exclusion.
