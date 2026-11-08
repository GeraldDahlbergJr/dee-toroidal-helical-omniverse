# Coupled Einstein–DEE evolution formulation

## Status

This document defines the next validation gate. It does **not** claim that the coupled evolution has passed.

## Evolution formulation

Use a BSSN-style 3+1 free-evolution formulation for the spacetime variables, coupled to the action-derived DEE matter fields `(rho, Theta, Psi)` and their normal derivatives `(Pi_rho, Pi_Theta, Pi_Psi)`.

The physical spatial metric and extrinsic curvature are decomposed as

`gamma_ij = exp(4 phi) gtilde_ij`,

`K_ij = exp(4 phi) (Atilde_ij + (1/3) gtilde_ij K)`.

The evolved gravitational variables are

- conformal factor `phi` (or an algebraically equivalent positive variable),
- conformal metric `gtilde_ij`, constrained by `det(gtilde)=1`,
- trace `K`,
- trace-free conformal extrinsic curvature `Atilde_ij`,
- conformal connection functions `Gamma_tilde^i`.

The matter source terms are constructed from the recovered nonlinear DEE stress tensor and projected as

`E = n^mu n^nu T_munu`,

`S_i = -gamma_i^mu n^nu T_munu`,

`S_ij = gamma_i^mu gamma_j^nu T_munu`,

`S = gamma^ij S_ij`.

No fixed-background substitution is allowed in this gate: the DEE fields source the gravitational evolution and the evolved geometry enters the DEE field equations.

## Gauge

Use standard moving-puncture-type gauge conditions as the initial implementation choice:

**1+log lapse**

`(partial_t - beta^i partial_i) alpha = -2 alpha K`.

**Gamma-driver shift**

`(partial_t - beta^j partial_j) beta^i = (3/4) B^i`,

`(partial_t - beta^j partial_j) B^i = (partial_t - beta^j partial_j) Gamma_tilde^i - eta B^i`.

`eta` is a documented numerical gauge parameter and must be included in every reproducibility record. Gauge changes require a separate checkpoint rather than silently replacing a failed run.

## Initial data

Start from the previously validated dynamical nonlinear initial-data checkpoint with nonzero matter momentum and nonzero extrinsic curvature. Convert its physical `(gamma_ij, K_ij)` into the BSSN variables without changing the physical initial data.

The initial lapse is `alpha=1` unless a separately documented pre-collapsed lapse is tested. The initial shift and driver are `beta^i=0`, `B^i=0` for the baseline run.

Before evolution begins, independently recompute the ADM Hamiltonian and momentum residuals from the reconstructed physical variables. The evolution is not started if this conversion destroys the previous constraint convergence.

## Independent diagnostics

At every output time independently evaluate

`H = R3 + K^2 - K_ij K^ij - 16 pi G E`,

`M^i = D_j(K^ij - gamma^ij K) - 8 pi G S^i`.

Also record:

- `L_inf` and RMS Hamiltonian residual,
- `L_inf` and RMS momentum residual,
- `min(alpha)`,
- `min(det(gamma_ij))`,
- conformal determinant error `|det(gtilde)-1|`,
- trace-free error `|gtilde^ij Atilde_ij|`,
- minimum and maximum `rho`,
- deviation of `rho` from the finite branch,
- representative curvature invariants where implemented,
- apparent/trapped-surface diagnostic when implemented,
- boundary-zone residuals separately from interior residuals.

## First coupled gate

Run at least three factor-of-two spatial resolutions over the same physical time interval and with the same physical initial data. Use a Courant factor that scales `dt` with grid spacing.

A provisional pass requires all of the following:

1. Matter and spacetime are both evolved; no fixed-background geometry is used.
2. The initial data retain the previously demonstrated Hamiltonian and momentum constraint convergence after conversion to evolution variables.
3. Hamiltonian and momentum residuals remain bounded over the stated interval and decrease with spatial refinement in the convergent regime.
4. The spatial metric remains positive definite and the reconstructed four-metric remains Lorentzian over the reported domain.
5. The finite DEE branch remains well-defined (`rho^2 > 0`) over the reported interval, or any loss is recorded as a physical/numerical failure rather than clipped away.
6. Algebraic BSSN constraints remain controlled.
7. Any trapped surface or curvature growth is reported, not filtered from the gate.

## Failure policy

The first configuration is retained in the record even if it fails. A later parameter or gauge scan is a separate experiment. Failure modes include constraint growth that does not improve under refinement, loss of metric regularity/signature, nonconvergent boundary contamination, or loss of the finite branch.

## Scope

Passing this gate would provide numerical evidence for survival of the **tested reduced coupled Einstein–DEE configuration over the tested time interval**. It would not by itself prove global nonlinear stability, a global toroidal solution, uniqueness, experimental validity, or the complete DEE framework.
