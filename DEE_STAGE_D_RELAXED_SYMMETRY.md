# Stage D: relaxed axial symmetry

Date: 2026-11-10 UTC. Exact parent: `514f5c77cb7250b064b54b3c989bde043c53ad70`.
Branch: `research/stage-d-relaxed-symmetry`. One agent; exact analysis only.

**Result — OBSTRUCTION IDENTIFIED under explicit global static hypotheses:**
removing axial translation invariance removes Stage C's pointwise constant-k
argument, but does not permit a regular isolated static soliton in the
nonnegative-potential, positive-kinetic sector described below. Two proofs and
their different hypotheses appear in the feasibility report. This is a broader
static no-go, not a theorem for every stationary or dynamical Einstein–DEE field.
No new existence result or complete spacetime is constructed.

## Baseline and unchanged covariant theory

Authority is `evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`,
Eqs. (3)–(9), (14)–(15). Read with Stage A½ symbolic verification, backreaction,
bounded response and boundary compatibility; both Stage B reports; both Stage C
reports and its verifier; and the matter existence analysis. Historical extended
matter existence, restricted localization verdicts, Eq. (34) correction, frozen
records and Run 43 retain their original scopes. The source action allows general
coupling functions; statements using the polynomial realization are marked.
Signature (-,+,+,+), G>0; smooth real fields and metric. Write

    S = integral sqrt(-g)[R/(16pi G) - (d rho)^2/2 - V
                                      - d_mu Q^T M d^mu Q/2],
    Q=(Theta,Psi)^T, M=[[K_Theta,Lambda],[Lambda,K_Psi]].

Variation before reduction, integration by parts with fixed boundary traces,
and symmetric inverse metric variation give

    E_rho = Box rho - V' - d_mu Q^T M' d^mu Q/2 = 0,
    E_Q = nabla_mu J^mu = 0, J^mu = M d^mu Q,
    T_mu nu = d_mu rho d_nu rho + d_mu Q^T M d_nu Q + g_mu nu L,
    G_mu nu = 8pi G T_mu nu.

The scalar boundary variation is
`-n_mu[d^mu rho delta rho + J^mu^T delta Q]`.
The off-shell identity is
`nabla_mu T^mu_nu = E_rho d_nu rho + E_Q^T d_nu Q`.
In particular no stress-conservation check alone replaces the amplitude EOM.
The polynomial realization, unchanged throughout, is

    M=I+rho² A, A=[[a,c],[c,b]],
    V=lambda_V(rho²-v²)²/4, lambda_V>0, v>0.

## Metric, symmetry and exact two-coordinate equations

Select a **static** subclass of stationary, orthogonally transitive axisymmetric
metrics, with no rotation or twist and symmetry-inheriting matter:

    ds²=-N(r,z)²dt²+H(r,z)²(dr²+dz²)+W(r,z)²dphi²,
    rho=rho(r,z), Q=Q(r,z), phi has period 2pi.

N,H>0 and W>0 off axis; a horizon is excluded in the global results. Isothermal
coordinates on the meridional metric are local; a global chart is an additional
assumption wherever coordinate boundary formulas are used. Invariant versions
of the proofs do not require one chart. W is not set to r: that would impose an
extra Einstein restriction. There is no z translation symmetry. Staticity is
an explicit restriction, not a consequence proved for arbitrary stationary
spacetimes. Non-inheriting matter with time-dependent phases and rotating metrics
is outside the selected ansatz.

Set F=NW, S=NH²W and sum a over r,z. The exact matter equations are

    E_rho = [partial_r(F rho_r)+partial_z(F rho_z)]/S - V'
             - [Q_r^T M' Q_r+Q_z^T M' Q_z]/(2H²)=0,
    S E_Theta = partial_r[F(K_Theta Theta_r+Lambda Psi_r)]
                 +partial_z[F(K_Theta Theta_z+Lambda Psi_z)]=0,
    S E_Psi = partial_r[F(Lambda Theta_r+K_Psi Psi_r)]
                 +partial_z[F(Lambda Theta_z+K_Psi Psi_z)]=0.

Equivalently multiply the last two by M inverse to obtain
`Q_rr+Q_zz+grad log F dot grad Q+M^{-1}M' (rho_r Q_r+rho_z Q_z)=0`.
The original amplitude equation contains `-rho sum Q_a^T A Q_a/H²` in the
polynomial case. These equations contain axial derivatives and amplitude–phase
couplings omitted in a radial reduction; they are not weak-field approximations.
For smooth phases `partial_z Q_r=partial_r Q_z`. If gradients are assigned as
independent unknowns this closure must accompany current divergence. In
particular `Q=k(r,z)z` has `Q_r=z k_r`, `Q_z=k+z k_z`:
there is no consistent prescription that changes only its axial derivative.

Locally current conservation permits axial flux to decrease while radial flux
increases: `partial_z(F M Q_z)=-partial_r(F M Q_r)`. It does not demand either
component be constant. The global integral obstruction below excludes nonzero
smooth phase gradients that taper to vacuum under its stated boundary conditions.

## Einstein equations and necessary constraints

Define h=H²(dr²+dz²)+W²dphi², D its connection, B=log H,
`C_ij=rho_i rho_j+Q_i^T M Q_j`,
`Kmat=h^{ij}C_ij`, `epsilon=Kmat/2+V`.
Static observers measure nonnegative energy for M>0,V>=0. Momentum density
j_i=0. With extrinsic curvature convention
`K_ij=(-partial_t h_ij+D_i beta_j+D_j beta_i)/(2N)`, beta=0 gives K_ij=0.
The ADM equations are therefore

    R[h]=16pi G epsilon,                 (Hamiltonian)
    D_j(K^{ij}-h^{ij}K)=8pi G j^i=0.     (momentum, identically)
    R[h]=-2/H² [B_rr+B_zz+(W_rr+W_zz)/W].

An arbitrary positive-energy profile on flat h fails the Hamiltonian constraint;
zero momentum is insufficient. For example a time-symmetric conformally flat
trial h=chi^4 delta must solve, with prescribed scalar maps,

    Delta_flat chi = -pi G chi C_flat -2pi G chi^5 V,
    C_flat=|grad_flat rho|²+sum grad_flat Q^T M grad_flat Q,
    chi>0, chi->1, regular axis/core.

This is a necessary elliptic problem, not a solved constraint construction.
Even solving it would produce initial data, not a static equilibrium. The
remaining static equations, using the trace-reversed four-dimensional Einstein
equations `R_mu nu=8pi G(C_mu nu+V g_mu nu)`, are

    D²N = -8pi G N V,
    R_ij[h] - N^{-1}D_iD_j N = 8pi G(C_ij+V h_ij).

They include the independent mixed radial–axial equation

    -W_rz/W+(B_z W_r+B_r W_z)/W
      -(N_rz-B_r N_z-B_z N_r)/N
        =8pi G(rho_r rho_z+Q_r^T M Q_z).

Here `D²N=[partial_r(W N_r)+partial_z(W N_z)]/(H²W)`.
The phi-phi equation, both meridional diagonal equations and lapse equation
must all hold; Bianchi relations do not license discarding unsolved equations.
The verifier independently computes the four-dimensional connection/Ricci
and checks the Hamiltonian geometry, lapse geometry, mixed equation and all
momentum components. Mixed scalar partials and both currents are checked too.

## Axis, vacuum, infinity and matching

A regular axis in this chart requires smooth even expansions in r for N,H,rho,Q,
W=r H(0,z)+O(r³), and N_r=H_r=rho_r=Q_r=0 at r=0.
The elementary-flatness condition is W_r/H=1 there. These imply no axis boundary
term for bounded phases and currents. A smooth filled axis forbids nonzero
circle-phase winding; M(v)>0 and M(0)=I do not allow an amplitude zero to hide
its singular energy. Scalar maps must be defined smoothly in regular Cartesian
charts, not merely in the degenerate cylindrical chart.

Strict positivity requires K_Theta>0 and det M>0 everywhere; in the polynomial
case det M=1+(a+b)rho²+(ab-c²)rho⁴. A positive semidefinite A suffices for
M>=I, but neither static theorem requires A positive. On any bounded amplitude
range, continuous strictly positive M has a uniform lower eigenvalue bound.
For unspecified functions nonnegative energy needs a separate V>=0 assumption.

In the positive sector, Stage C vacuum rigidity still applies without translation
symmetry: T=0 implies constant rho,Q, V=V'=0. Polynomial vacuum is rho=v
(or -v for a signed scalar), constant phases; rho=0 has positive potential.
Necessary isolated-object requirements are finite `E_m=integral epsilon dvol_h`,
T in an orthonormal frame tending to zero, rho tending to a permitted vacuum,
phase gradients tending to zero, no flux sources, and appropriate gravitational
falloff. Finite integrated energy alone does not imply pointwise decay without
additional regularity/falloff assumptions. E_m is proper matter energy, not ADM
mass; redshift-weighted Killing energy is `integral N epsilon dvol_h`.

For the global results require a complete smooth slice with one asymptotically
Euclidean end, no inner boundaries or singularities, N>0, N->1,
h_ij=delta_ij+O(1/R), partial h=O(R^-2), N=1-G M_K/R+o(R^-1),
and the differentiated falloff/integrability needed for the ADM and Komar limits
and their equality in a static Einstein solution. R is asymptotic spherical
radius sqrt(r²+z²), not cylindrical r at fixed z. Standard isolated falloff may
be chosen `rho-v=O(R^-p)`, `Q-Q_infty=O(R^-q)` with derivatives one power
faster, p>3/2 and q>1/2; bounded smooth M near vacuum then suffices for matter
energy integrability and vanishing integral-identity boundary terms. These
are sufficient power bounds, not claimed necessary optimal rates or solutions.
A nonzero conserved axial k is incompatible with these conditions.

Exact finite vacuum matching instead requires continuous scalar traces and
normal scalar derivative/phase currents, plus equal induced metric and extrinsic
curvature (no shell). Smooth joins are stronger. No radial ODE uniqueness
argument is imported as a PDE proof. Under the global theorems in the companion
report, smooth compact support is excluded as a special case; on bounded or
multiply connected domains without those hypotheses, a separate elliptic
unique-continuation argument would be needed. Asymptotic tails have no finite
vacuum interface and must be tested against the global identities separately.
