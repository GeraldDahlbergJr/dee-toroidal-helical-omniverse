# Coupled dynamical nonlinear DEE constraint gate

## Purpose

This gate removes the time-symmetric shortcut used by the preceding conformal
initial-data checkpoint.  It must not be reported as passed until both the
Hamiltonian and momentum equations are solved with nonzero matter momentum.

## Data and decomposition

Use the recovered DEE fields `(rho, Theta, Psi)` and canonical normal
velocities `(Pi_rho, Pi_Theta, Pi_Psi)` to construct the Eulerian energy and
momentum densities `E` and `S_i` from the action-derived stress tensor.
Require at least one nonzero `Pi` and a spatial phase gradient so that
`max |S_i| > 0`.

Adopt the conformal transverse-traceless decomposition

`gamma_ij = psi^4 delta_ij`, `K = 0`,

`K^ij = psi^-10 (L W)^ij`,

where

`(L W)_ij = d_i W_j + d_j W_i - (2/3) delta_ij d_k W^k`.

The unknowns are therefore the conformal factor `psi` and vector potential
`W^i`; neither the physical spatial metric nor extrinsic curvature is
prescribed after the matter data are selected.

## Coupled constraints

With the stated flat conformal background and maximal slicing, solve the
Lichnerowicz-York system

`Delta psi + (1/8) psi^-7 Abar_ij Abar^ij + 2 pi G psi^5 E = 0`,

`Delta_L W^i = 8 pi G psi^6 S^i`,

using a single, documented convention for the conformal scaling of matter
sources.  The physical solution must then be reconstructed and checked using
the ADM equations directly, rather than judging success only from the
elliptic equations used by the solver.

## Independent residuals

After solving, independently recompute

`H = R3 + K^2 - K_ij K^ij - 16 pi G E`,

`M^i = D_j (K^ij - gamma^ij K) - 8 pi G S^i`.

The residual evaluator must be separate from the nonlinear solve path to
avoid an implementation cancelling its own discretization error.

## Refinement gate

Run at least four factor-of-two refinements.  Record for both constraints:

- maximum absolute residual,
- RMS residual,
- observed convergence order,
- boundary residuals,
- `min(psi)` and metric signature/positivity diagnostics,
- `max |S_i|` to prove the gate is genuinely dynamical,
- `max |K_ij|` to prove the extrinsic-curvature sector is active.

A pass requires:

1. `max |S_i| > 0` and `max |K_ij| > 0`;
2. positive regular conformal factor throughout the domain;
3. Hamiltonian and momentum residuals decrease systematically with refinement;
4. the asymptotic observed order is compatible with the documented finite-
   difference order (target 2 for the current numerical stack);
5. no constraint is declared solved merely because a symmetry choice makes
   its source vanish.

## Failure policy

If the coupled elliptic solve has no regular solution for the chosen DEE data,
or either independently evaluated constraint fails to converge toward zero,
that data set fails this gate.  Parameters may subsequently be scanned, but
the failed case must remain recorded rather than silently replaced.

## Scope

Passing this gate establishes constraint-satisfying *initial data* for the
specified reduced conformal sector.  It does not establish nonlinear
stability or survival under time evolution.  Evolution with constraint,
curvature, null-cone, signature, trapped-surface, and DEE separation
monitoring is a separate subsequent gate.
