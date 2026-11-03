# Nonlinear DEE source for 3+1 Einstein evolution

With signature (-,+,+,+), let `Pi_A = -n^mu partial_mu phi_A`,
`D_i phi_A = partial_i phi_A`, and let the field-space kinetic matrix be

`F_AB = [[1,0,0],[0,1+a*rho^2,c*rho^2],[0,c*rho^2,1+b*rho^2]]`

for `phi_A=(rho,Theta,Psi)`. The matter Lagrangian is
`L = (F_AB Pi_A Pi_B - F_AB gamma^ij D_i phi_A D_j phi_B)/2 - V`, where
`V=lambda_V (rho^2-v^2)^2/4`. Direct metric variation gives
`T_mu nu = F_AB partial_mu phi_A partial_nu phi_B + g_mu nu L`.
The full projections are therefore

`E = (F_AB Pi_A Pi_B + F_AB gamma^ij D_i phi_A D_j phi_B)/2 + V`,

`S_i = F_AB Pi_A D_i phi_B`,

`S_ij = F_AB D_i phi_A D_j phi_B + gamma_ij L`, and
`S = gamma^ij S_ij`.

Repeated field indices are summed. In particular, the mixed terms are
`c*rho^2*(Pi_Theta Pi_Psi + DTheta.DPsi)` in E,
`c*rho^2*(Pi_Theta D_i Psi + Pi_Psi D_i Theta)` in S_i, and
`c*rho^2*(D_i Theta D_j Psi + D_i Psi D_j Theta)` in the dyadic part of
S_ij. Its pressure term includes `c*rho^2*(Pi_Theta Pi_Psi - DTheta.DPsi)`.
The lapse and shift are already accounted for in `Pi_A`.

These are the source inputs for a BSSN implementation. In standard ADM
conventions the matter part of the K_ij RHS is
`-8*pi*G*alpha*(S_ij - gamma_ij*(S-E)/2)`; a conformal trace-free
BSSN RHS uses its trace-free part, while the K trace RHS contains
`4*pi*G*alpha*(E+S)` in the usual constraint-substituted form. Gauge,
curvature derivatives, boundary treatment, and geometric RHS still need
implementation and convergence verification. No coupled evolution is
claimed by this source projection.

`tests/test_matter_projections.py` independently reconstructs the full
four-metric and covariant tensor for nontrivial lapse, shift, off-diagonal
spatial metric, and nonzero mixed coupling, then contracts with the normal.
