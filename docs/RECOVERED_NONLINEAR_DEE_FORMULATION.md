# Recovered nonlinear DEE formulation

## Source status

This record transcribes the nonlinear matter formulation and nonlinear initial-value target from the author's archived DEE manuscript, rather than inventing a new nonlinear theory. The separate helical-perturbation draft supplies controlled perturbation variables and diagnostics but explicitly leaves the final background metric/matter model to the DEE formulation.

## Conventions

Coordinates: `x^mu = (t,r,phi,z)`. Signature: `(-,+,+,+)`. Units: `c=1`.

Fields:

- `rho(x)`: finite-separation / collapse scalar.
- `Theta(x)`: toroidal transport phase.
- `Psi(x)`: helical transport phase.

## Covariant action

S = integral d^4x sqrt(-g) [ R/(16 pi G)
  - 1/2 (grad rho)^2 - V(rho)
  - 1/2 K_Theta(rho) (grad Theta)^2
  - 1/2 K_Psi(rho) (grad Psi)^2
  - Lambda(rho) grad_mu Theta grad^mu Psi ].

This action is the model ansatz. Its consequences are to be tested; its microscopic origin is not assumed.

## Euler-Lagrange equations

Radial scalar:

box rho - V'(rho)
- 1/2 K_Theta'(rho) (grad Theta)^2
- 1/2 K_Psi'(rho) (grad Psi)^2
- Lambda'(rho) grad_mu Theta grad^mu Psi = 0.

Transport currents:

nabla_mu [ K_Theta grad^mu Theta + Lambda grad^mu Psi ] = 0,

nabla_mu [ K_Psi grad^mu Psi + Lambda grad^mu Theta ] = 0.

## Stress-energy

T_mn = grad_m rho grad_n rho
+ K_Theta grad_m Theta grad_n Theta
+ K_Psi grad_m Psi grad_n Psi
+ Lambda (grad_m Theta grad_n Psi + grad_m Psi grad_n Theta)
+ g_mn L_m,

with

L_m = -1/2(grad rho)^2 - V
      -1/2 K_Theta(grad Theta)^2
      -1/2 K_Psi(grad Psi)^2
      -Lambda grad_a Theta grad^a Psi.

The nonlinear gravitational equation is

G_mn[g] = 8 pi G T_mn[g,rho,Theta,Psi].

## Minimal polynomial realization

V(rho) = lambda_V/4 (rho^2-v^2)^2,
K_Theta = 1 + a rho^2,
K_Psi = 1 + b rho^2,
Lambda = c rho^2.

Define

X_Theta = (grad Theta)^2,
X_Psi = (grad Psi)^2,
X_ThetaPsi = grad_mu Theta grad^mu Psi.

The finite branch satisfies

rho_* [lambda_V(rho_*^2-v^2) + a X_Theta + b X_Psi + 2c X_ThetaPsi] = 0,

hence

rho_*^2 = v^2 - (a X_Theta+b X_Psi+2c X_ThetaPsi)/lambda_V.

For a positive branch the homogeneous radial curvature is

m_eff^2 = 2 lambda_V rho_*^2 > 0.

The transport kinetic matrix is

K = [[K_Theta,Lambda],[Lambda,K_Psi]],

and positive definiteness requires `K_Theta>0` and `det(K)>0` (with `K_Psi>0` following for a symmetric 2x2 positive-definite matrix). For every null k,

T_mn k^m k^n = (k.grad rho)^2 + u^T K u,

where `u=(k.grad Theta,k.grad Psi)`. Therefore positive semidefinite K is sufficient for the NEC of this matter action.

## Explicit local transport profiles

Theta = omega_Theta t + m_Theta phi + k_Theta z,
Psi   = omega_Psi t   + m_Psi phi   + k_Psi z.

These are local non-oscillatory gradient profiles. Global angular single-valuedness remains a separate regularity/winding condition.

## Helical perturbation target

The archived helical draft defines

Phi(z,phi,t) = k_z z + m phi - omega t + Phi_0,

delta X = A F(rho) W(Phi),

with an initially sinusoidal `W(Phi)|t=0 = cos(Phi)` but nonlinear evolution permitted to generate higher harmonics. Its parameter set is `{A,k_z,m,omega,chi}`. These are perturbation/initial-data controls, not replacements for the covariant DEE action above.

## Nonlinear initial-value target from the DEE manuscript

The weak-field solution is explicitly proposed as a seed:

g_mn(t0,r) = eta_mn + h_mn(r),
partial_t g_mn(t0,r) = dot(A_DEE)(t0) h_mn^(DEE)(r),
rho(t0,r) = rho_* + delta rho(r),
partial_t rho(t0,r) = delta dot(rho)(r),

with corresponding Theta/Psi data.

These seed data are **not automatically constraint satisfying**. They must be adjusted/solved so that

(3)R + K^2 - K_ij K^ij = 16 pi G E,
D_j(K^ij-gamma^ij K) = 8 pi G S^i.

## Survival criteria

A nonlinear run can count as a survival result only over a stated simulation interval and resolution hierarchy if:

1. `rho_min(t)>0`;
2. the full-null-cone minimum `T_mn k^m k^n >= 0` if NEC preservation is part of the tested claim;
3. curvature invariants remain finite;
4. Hamiltonian and momentum constraint residuals converge toward zero under refinement;
5. metric signature and regularity remain admissible;
6. trapped-surface/horizon diagnostics are explicitly monitored where applicable;
7. the weak-field limit recovers the independently verified linearized sector.

Failure of any required criterion is to be reported as failure of the tested parameter set, not tuned away after inspection.

## Important boundary

The archived helical perturbation draft does not itself specify a final nonlinear spacetime metric ansatz. It explicitly identifies its parameter map as initial-data targeting for a later Einstein-matter calculation. The nonlinear covariant content available for validation is therefore the DEE Einstein + three-scalar action above, together with constraint-satisfying initial data to be constructed from the weak-field seed.