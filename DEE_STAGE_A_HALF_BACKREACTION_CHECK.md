# Stage A½: matter backreaction and conservation check

Date: 2026-11-10 UTC. Base: `bd32e856b49bc0c227d235ef68f345d72acb2b21`.
Conventions, recovered inputs and action sources are recorded in `DEE_STAGE_A_HALF_SYMBOLIC_VERIFICATION.md`. All results below are local first-order statements in the original action.

## Variation of the original action

Write M'=2rho A. Integration by parts in the scalar variations and inverse-metric variation give

    E_rho = Box_g rho - V'(rho)
            -rho g^{mu nu}(partial_mu Q)^T A partial_nu Q,
    E_Q = nabla_mu J^mu,
    J^mu = M(rho) partial^mu Q,
    T_mu nu = partial_mu rho partial_nu rho
              +(partial_mu Q)^T M(rho) partial_nu Q + g_mu nu L_m.

The two phase equations, explicitly, are

    nabla_mu[(1+a rho²) partial^mu Theta+c rho² partial^mu Psi]=0,
    nabla_mu[c rho² partial^mu Theta+(1+b rho²) partial^mu Psi]=0.

Here E=0 defines each Euler–Lagrange equation; these signs fix the residual convention. Symmetry of A and commutation of scalar second derivatives imply the off-shell identity

    nabla_mu T^mu_nu = E_rho partial_nu rho + (partial_nu Q)^T E_Q.

For example differentiating the kinetic stress gives its Euler–Lagrange divergence plus derivatives of the kinetic energy; the latter cancel against partial_nu L_m, leaving exactly the expression above. Thus on-shell covariant conservation follows from the three scalar equations.

## Linearized amplitude

Introduce a formal parameter tau with g=eta+tau h, rho=rho_0+tau u and Q=kz+Q_0+tau f. Retain the coefficient of tau, with u=delta rho and f=(delta Theta,delta Psi)^T. Then

    g^{mu nu}=eta^{mu nu}-tau h^{mu nu}+O(tau²),
    sqrt(-g)=1+tau h/2+O(tau²),
    delta[(partial Q)^T A (partial Q)]
        =-C h_zz+2 k^T A partial_z f.

The last expression denotes the contracted spacetime kinetic scalar. Since rho_0 is constant, variation of Box_g acting on the background contributes zero. Also V''(rho_0)+C=2lambda_V rho_0² by the background branch relation. Hence

    E_rho^(1) = Box u - 2lambda_V rho_0² u
                +rho_0 C h_zz -2rho_0 k^T A partial_z f.

**VERIFIED:** the supplied amplitude equation has the correct sign and coefficient:

    Box delta rho - 2lambda_V rho_0² delta rho
       = -rho_0 C h_zz +2rho_0 k^T A partial_z delta Q.

The term with A is a row-vector contraction, not an additional scalar coupling. This derivation does not require harmonic gauge. The positive mass squared concerns the amplitude with prescribed phase gradients; it does not prove stability of every coupled mode.

For fixed phases f=0 and a static, z-independent u(x,y), Box u=Delta_perp u. Therefore the supplied transverse equation is **VERIFIED**:

    (Delta_perp-2lambda_V rho_0²)u
       =2 pi G rho_0 C(B+2epsilon)s.

No boundary data or selected response solution is supplied or solved here. Admissible matter backreaction with global boundaries remains **NOT YET ESTABLISHED**.

## Frozen-matter residual and its sign

Set u=f=0. Using ordinary h, not bar h, the original amplitude EOM has residual

    E_rho^(1)|frozen = rho_0 C h_zz
       = -2 pi G rho_0 C(B+2epsilon)s.

**VERIFIED:** the previously reported residual is correct under this EOM convention. Equivalently it is the reduced amplitude LHS minus RHS. The forcing RHS itself is the opposite sign, +2 pi G rho_0 C(B+2epsilon)s. Defining residual as RHS minus LHS would reverse the reported sign.

**CORRECTION REQUIRED:** frozen matter is generically inconsistent with its amplitude equation in the candidate metric. The obstruction vanishes identically if G rho_0 C(B+2epsilon)=0; it vanishes pointwise on s=0 without vanishing on a neighborhood. With G>0, rho_0>0, positive M, lambda_V>0 and nonzero k, B+2epsilon=2B+2V>0, leaving C=0 as the coefficient exception. No universal sign for C is assumed for arbitrary symmetric A. Positive-definite A and nonzero k imply C>0 and a negative residual for s>0. None of these exceptions constitutes a global spacetime construction.

## Both linearized phase-current equations

Put M_0=M(rho_0), h=eta^{mu nu}h_mu nu and raise metric indices with eta. The variation of the current density is

    delta(sqrt(-g) J^mu)
      =M_0 partial^mu f+2rho_0 u A k delta^mu_z
       +(h/2 eta^{mu z}-h^{mu z})M_0 k.

Define D^nu=partial_mu h^{mu nu}-(1/2)partial^nu h. Then

    E_Q^(1)=M_0 Box f+2rho_0 A k partial_z u-M_0 k D^z.

The verifier differentiates the current density for two arbitrary phase perturbations and an arbitrary symmetric spacetime-dependent h; both residuals relative to this formula are exactly zero. In harmonic gauge D^nu=0 the equations become

    (1+a rho_0²) Box delta Theta+c rho_0² Box delta Psi
      +2rho_0(a k_Theta+c k_Psi)partial_z u=0,
    c rho_0² Box delta Theta+(1+b rho_0²) Box delta Psi
      +2rho_0(c k_Theta+b k_Psi)partial_z u=0.

**VERIFIED:** fixed phases and transverse u satisfy both at this order in the supplied metric. This does not solve the amplitude response or its boundaries.

## First-order stress and four conservation identities

Let d=k^T M_0 partial_z f. Background stationarity makes the coefficient of u in the variation of V+(1/2)(partial Q)^T M(partial Q) cancel: V'(rho_0)+rho_0 C=0. Consequently

    delta T_mu nu = 2rho_0 C u delta_mu_z delta_nu_z
       +delta_mu_z k^T M_0 partial_nu f
       +delta_nu_z k^T M_0 partial_mu f
       -epsilon h_mu nu-eta_mu nu[d-(B/2)h_zz].

Terms quadratic in amplitude gradients are second order. The symbolic verifier independently extracts all entries of delta T from the covariant action stress and checks this expression exactly.

For constant T^(0), define

    delta T^mu_nu = eta^{mu alpha}delta T_alpha nu
                   -h^{mu alpha}T_alpha nu^(0),
    Gamma^(1)lambda_mu nu = (1/2)eta^{lambda alpha}
       (partial_mu h_alpha nu+partial_nu h_alpha mu-partial_alpha h_mu nu),
    C_nu = partial_mu delta T^mu_nu
           +Gamma^(1)mu_mu lambda T^(0)lambda_nu
           -Gamma^(1)lambda_mu nu T^(0)mu_lambda.

The background is on-shell and partial_nu rho_0=0, so the action identity linearizes to

    C_t=0, C_x=0, C_y=0,
    C_z=k^T E_Q^(1)
       = k^T M_0 Box f+2rho_0 C partial_z u-B D^z.

**VERIFIED:** all four identities hold off-shell in u,f,h at first order, using only the background relation. The connection and index-raising terms are essential; partial^mu delta T_mu nu alone is not the covariant conservation residual. The verifier builds the Christoffel terms with arbitrary h and checks each identity separately.

For the candidate h with frozen u=f=0, both phase residuals and all four covariant conservation residuals are exactly zero, while E_rho^(1) generally is not. Amplitude EOM failure first enters the off-shell conservation identity through E_rho^(1)partial_nu u at second order, since the background amplitude gradient vanishes. **CORRECTION REQUIRED:** first-order stress conservation cannot certify amplitude backreaction in this background.

## Perturbative meaning and preservation

The formal tau expansion tests matter response to a weak prescribed metric. For gravitational coupling counting with fixed background stress, h,u,f=O(G), a regular response gives delta T=O(G) and G delta T=O(G²). Thus sourcing G^(1)[h] with T^(0) is the leading weak-gravity equation; it does not assert that frozen fields solve their own O(G) equations or that this is an exact linearization about a flat solution of the full coupled system with nonzero T^(0).

Reproduce with `python DEE_STAGE_A_HALF_VERIFIER.py` (Python 3; checked with SymPy 1.14.0). All exact comparison residuals print zero, followed by the generally nonzero frozen amplitude expression and the opposite-sign transverse forcing. No numerical solving is involved.

Standing **MATTER EXISTENCE: ESTABLISHED**, **LOCALIZATION: OBSTRUCTION IDENTIFIED** within the original stated scope, and **COMPLETE EINSTEIN–DEE SPACETIME: NOT YET ESTABLISHED** are preserved. Existing action, Eq. (34) correction, analytical evidence, manuscripts and Run 43 artifacts remain unchanged. Stop after publication and remote verification of this evidence package; no nonlinear construction or numerical solving follows.
