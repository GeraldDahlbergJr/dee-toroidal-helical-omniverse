# Stage E: stationary Einstein–DEE equations

Date: 2026-11-10 UTC. Exact Stage D parent:
`570a9fdcf8ce7fb62e8d377bb535e6da6b42595d`.
Branch: `research/stage-e-stationary-analysis`. One Codex agent, exact analysis;
no numerical PDE, evolution, parameter search, manuscript or Stage F.

**VERIFIED:** the equations below follow from the original action. **CONDITIONAL:**
the global hypotheses in the companion no-go report are additional assumptions,
not consequences of stationarity. **OBSTRUCTION IDENTIFIED:** under those
hypotheses the selected stationary class reduces to the Stage D vacuum result.
**NOT YET ESTABLISHED:** an obstruction for every stationary stress realization
or every stationary spacetime of the original theory.

## Preserved inputs and conventions

Read the original source
`evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`, Eqs. (3)–(9),
(14)–(15), Stage A½ symbolic, backreaction, bounded-response and boundary reports,
both Stage B reports, both Stage C reports, and both Stage D reports. Earlier
extended-matter existence, perturbative collar, localization and provenance
verdicts keep their stated scopes; the Eq. (34) correction remains required.
All pre-existing files, frozen verification evidence and Run 43 are preserved.

Signature (-,+,+,+), G>0, c=1. The unchanged theory is

    S = integral sqrt(-g) [R/(16 pi G) + L_m],
    L_m = -(d rho)^2/2 - V(rho) - d_mu Q^T M(rho) d^mu Q/2,
    Q=(Theta,Psi)^T, A=[[a,c],[c,b]], M=I+rho² A,
    V=lambda_V (rho²-v²)²/4, lambda_V>0, v>0.

A is real symmetric; no sign for A or M' is inferred from M>0.
Require M(rho)>0 on the entire solution. Equivalently 1+a rho²>0 and
`1+(a+b)rho²+(ab-c²)rho⁴>0`. M(0)=I; M(v)>0 is essential below.
The nonnegative separation interpretation of rho is retained when desired;
the gravitational obstruction also allows a signed scalar approaching +v.

Variation before symmetry reduction, with fixed boundary traces, gives

    E_rho = Box rho - V' - d_mu Q^T M' d^mu Q/2 = 0,
    E_Q = nabla_mu J^mu = 0, J^mu=M d^mu Q,
    T_mu nu = C_mu nu + g_mu nu L_m,
    C_mu nu = d_mu rho d_nu rho + d_mu Q^T M d_nu Q,
    G_mu nu=8 pi G T_mu nu,
    R_mu nu=8 pi G(C_mu nu+V g_mu nu).

The scalar surface variation is
`-n_mu[d^mu rho delta rho + (J^mu)^T delta Q]`.
The exact off-shell identity is
`nabla_mu T^mu_nu=E_rho d_nu rho+E_Q^T d_nu Q`.
This identity does not substitute for any scalar equation.

## Selected symmetry class, inverse metric and matter equations

Take commuting stationary and axial Killing fields k=partial_t,
eta=partial_phi, phi of period 2pi, and **assume circularity** (orthogonal
transitivity). In adapted charts, with meridional metric gamma_AB positive,

    ds²=-N²dt²+gamma_AB dx^A dx^B+W²(dphi-omega dt)²,
    rho=rho(x^A), Q=w t+m phi+q(x^A),
    U=w+omega m, sigma=sqrt(det gamma), F=NW.

w,m are constant two-vectors; the phase gradients, not necessarily the phases,
inherit both Killing symmetries. The stress inherits them. This is a selected
field representation of stationary stress, not an exhaustive classification.
N>0 and W>0 off axis; g_tt may become positive, so an ergoregion has not been
excluded at the level of the local equations. A global slice and time coordinate
are separate hypotheses for the theorems.

    g^tt=-1/N², g^tphi=-omega/N²,
    g^phiphi=1/W²-omega²/N², g^AB=gamma^AB,
    det g=-N² W² det gamma, sqrt(-g)=N W sigma.

For any symmetric field-space matrix B define

    K[B]=gamma^AB q_A^T B q_B + m^T B m/W² - U^T B U/N².

The exact amplitude equation is

    [partial_A(NW sigma gamma^AB rho_B)]/(NW sigma)
        -lambda_V rho(rho²-v²) -rho K[A] = 0.

Here K[A] uses the coupling matrix A; meridional indices are contracted inside
K[B]. Equivalently write its last term as `-K[M']/2`.
The sign of the temporal contribution is positive in E_rho. No additional
centrifugal term is inserted. The two-component currents and equations are

    J^t=-M U/N²,
    J^phi=M m/W²-omega M U/N²,
    J^A=M gamma^AB q_B,
    partial_A(NW sigma M gamma^AB q_B)=0.

The first component uses `(1+a rho²)Theta_B+c rho² Psi_B`; the second uses
`c rho² Theta_B+(1+b rho²)Psi_B`. Both equations must hold. Temporal and
angular current divergences vanish because their coefficients are independent
of t,phi, not because these currents vanish. Mixed scalar derivatives commute;
assigning gradients independently also requires closure d(dQ)=0.

## Complete stress and circularity compatibility

Define

    P=U^T M U/N²,
    X=gamma^AB(rho_A rho_B+q_A^T M q_B)+m^T M m/W²,
    L_m=(P-X)/2-V.

Every independent covariant component is

    T_tt=w^T M w+(-N²+W² omega²)L_m,
    T_tphi=w^T M m-W² omega L_m,
    T_phiphi=m^T M m+W² L_m,
    T_tA=w^T M q_A,
    T_phiA=m^T M q_A,
    T_AB=rho_A rho_B+q_A^T M q_B+gamma_AB L_m.

The unit slice normal is n=(partial_t+omega partial_phi)/N. For
h=gamma_AB dx^A dx^B+W²dphi² and beta=-omega partial_phi,

    epsilon=T(n,n)=(P+X)/2+V,
    j_A=-T(n,partial_A)=-U^T M q_A/N,
    j_phi=-T(n,partial_phi)=-U^T M m/N,
    S_ij=C_ij+h_ij L_m,
    S=trace_h S=(3P-X)/2-3V,
    epsilon+S=2P-2V.

The orthonormal angular derivative is m/W; the temporal derivative is U/N.
Positive kinetic contributions cannot cancel one another in epsilon. In
particular T_hat0hatphi=U^T M m/(NW), whereas j_hatphi is its negative
under the ADM convention above. The coordinate T_tphi can be nonzero even
if j_phi=0 because of the metric-times-L term.

For this block metric, R_tA=R_phiA=0 identically. Einstein equations therefore
**require**

    w^T M q_A=0, m^T M q_A=0 (each A),

or equivalently U^T M q_A=0 and m^T M q_A=0. These are actual mixed field
equations, not optional matter restrictions or automatic current identities.
If w,m span the two-dimensional phase space they force q_A=0, since M is
invertible; if dependent they need not. Circularity itself was imposed.
A noncircular solution could instead have nonzero meridional stress fluxes.
No circularity theorem for arbitrary DEE stationary stress is asserted.

## Einstein equations and ADM constraints, including the shift

Use the same convention as Stage D:

    K_ij=(-partial_t h_ij+D_i beta_j+D_j beta_i)/(2N).

D is the h connection, indices i,j run over (x^1,x^2,phi). Direct computation:

    K_AB=K_phiphi=0, K_Aphi=-W² partial_A omega/(2N), K=0,
    K_ij K^ij = W² gamma^AB omega_A omega_B/(2N²) =: K2.

All constraints are

    R[h]-K2=16 pi G epsilon,                        (Hamiltonian)
    D_j K^j_i=8 pi G j_i.                           (momentum)

The meridional left sides vanish and enforce j_A=0. For reference, if
nabla_gamma and Delta_gamma are the two-dimensional connection and Laplacian,

    R_AB[h]=R_AB[gamma]-W^-1 nabla_gamma_A nabla_gamma_B W,
    R_phiphi[h]=-W Delta_gamma W, R_Aphi[h]=0,
    R[h]=R[gamma]-2 W^-1 Delta_gamma W,
    D²N=(W sigma)^-1 partial_A(W sigma gamma^AB N_B).

These formulas retain an arbitrary meridional metric. In local isothermal
coordinates gamma=H²(dr²+dz²), one may substitute
`R[h]=-2/H²[(log H)_rr+(log H)_zz+(W_rr+W_zz)/W]`;
W=r is not an available arbitrary gauge choice.

The trace evolution equation, including K2, is

    D²N=N K2+4 pi G N(epsilon+S)
        =N K2+8 pi G(U^T M U/N-N V).

All six independent spatial Einstein equations are equivalently

    0=-D_i D_j N+N(R_ij[h]-2 K_ik K^k_j)
          +(Lie_beta K)_ij-8 pi G N(C_ij+V h_ij).

Together with the constraints and stationary definition of K these are a
complete ADM representation of Einstein's equations; the displayed lapse
trace is redundant when that full system holds. In particular include both
meridional diagonal equations, their off-diagonal equation, phi-phi and A-phi.
For unambiguous component use

    (Lie_beta K)_ij=beta^k partial_k K_ij
                       +K_kj partial_i beta^k+K_ik partial_j beta^k,
    (Lie_beta K)_AB=W² omega_A omega_B/N,
    (Lie_beta K)_Aphi=(Lie_beta K)_phiphi=0.

For example the meridional equations reduce to

    R_AB[h]-N^-1 D_A D_B N+W² omega_A omega_B/(2N²)
        =8 pi G(rho_A rho_B+q_A^T M q_B+V gamma_AB),

and the angular equation is

    R_phiphi[h]-N^-1 D_phi D_phi N
            -W^4 gamma^AB omega_A omega_B/(2N²)
        =8 pi G(m^T M m+V W²).

The A-phi equations impose m^T M q_A=0. Invariantly
`R_nn=D²N/N-K2=8 pi G(P-V)` and
`R_ni=-D_j K^j_i= -8 pi G j_i` fix the momentum sign. This explicitly
accounts for stationary shift and extrinsic curvature; the Stage D static
lapse equation cannot be used while K2 or P is nonzero.

## Frame dragging and angular momentum

Let `a^{AB}=sigma W³ gamma^AB/N`. The phi momentum equation is

    D_j K^j_phi=-(2 W sigma)^-1 partial_A(a^{AB} omega_B)
               =8 pi G j_phi,
    partial_A(a^{AB} omega_B)
          =16 pi G sigma W U^T M m/N.

This is the frame-dragging equation with the conventions above. It is coupled
to the other Einstein and matter equations, and is not a prescribed rotation
profile. Multiplication by omega yields the exact local identity

    partial_A(omega a^{AB} omega_B)
      =a^{AB} omega_A omega_B
         +16 pi G sigma W omega U^T M m/N.

No sign for the last term is available for arbitrary w,m. If w=m=0, j_phi=0
and the last term vanishes. Staticity then depends on the full boundary
identity proved in Theorem E3, not just on the pointwise source vanishing.

With eta=partial_phi and outward normal s on a large spatial sphere,

    J_ADM=(1/(8 pi G)) lim integral K_ij eta^i s^j dA
         =-(1/(16 pi G)) lim integral W² N^-1 D_s omega dA.

The stationary angular Komar charge agrees with this under standard falloff.
Contracting the momentum constraint with the spatial Killing field eta gives
`D_j(K^j_i eta^i)=8 pi G j_i eta^i`. Hence J_ADM equals the volume integral
of j_phi only if all other boundary fluxes vanish and the volume limit exists.
For omega=2GJ/R³+o(R^-3) this surface formula has the stated positive J sign.
An inner horizon can carry angular momentum even in matter vacuum.

## Axis and sufficient asymptotic conditions

**CONDITIONAL:** a smooth filled axis is part of the regular sector studied.
In proper distance s from it, W=s+O(s³) (in an adapted elementary-flat chart),
N and omega are bounded smooth rotational scalars with even expansions.
Smooth invariant rho,q have vanishing transverse first derivative there.
The degeneracy W=0 is a coordinate degeneracy, not an excised physical boundary.
M is continuous and strictly positive at each regular axis point.

For the global result choose a complete orientable smooth slice with one
asymptotically Euclidean end, no horizons, singularities, excisions or other
inner boundaries, a global circular decomposition and N>0. On the end take
h_ij=delta_ij+O(R^-1), differentiated standard ADM falloff, N=1+O(R^-1),
DN=O(R^-2), omega=O(R^-3), Domega=O(R^-4), and corresponding regular angular
bounds. Here R is spherical distance and W~R sin(theta), not fixed-r z-infinity.
These are sufficient isolated rotation falloffs, not conclusions of the local
PDE. Require finite proper matter energy, suitable integrability of the bulk
identities, vacuum-frame stress decay, rho->v and M->M(v)>0. Once w=m=0,
sufficient scalar rates are rho-v=O(R^-p), D rho=O(R^-p-1), p>3/2,
and local phase derivatives O(R^-q-1), q>1/2. Where phase lifts are used,
q-q_infty=O(R^-q) suffices. The gravitational theorem needs no global lifts.
Assume also the standard differentiated stationary falloff ensuring Komar/ADM
mass equality; these are the same mass hypotheses as Stage D after staticity.

## Machine verification scope

`python DEE_STAGE_E_VERIFIER.py` independently reduces the action for all three
scalar equations and varies the inverse metric for all ten stress components.
It checks both currents, all four off-shell stress divergences with an explicit
connection, the inverse/determinant, observer energy and angular momentum,
active mass, direct four-dimensional Ricci versus independent three-dimensional
ADM geometry, every shift-defined K component, Hamiltonian, three momentum
identities, lapse identity, all six spatial curvature identities, mixed
circular components, rotation PDE and Green product rule, and asymptotic and
winding algebra. The local chart uses arbitrary functions H,N,W,omega,rho,q
of two variables. Local isothermal coordinates are legitimate in two dimensions;
no global coordinate coverage is inferred. Inequalities, boundary limits,
completeness, topology and the positive mass theorem are mathematical arguments
in the companion report, not outputs certified by SymPy.
