# Run 43: W7-X-inspired toroidal magnetic potential pilot

**Status:** New independent pilot; no modification to frozen Run 43 baseline or its 12/12 audit.

## Goal
Determine whether flux-compatible magnetic-potential variables and commuting discrete derivatives can provide a consistent interface to the existing DEE toroidal/helical perturbation model. This is **not yet an explanation or repair** for the experimental inner-momentum 9/12 failures.

## Geometry and equations
Use coordinates (r, theta, phi) on a regular toroidal chart and the one-form
A = 0.3 r^2 dtheta + [-0.2 r^2 + epsilon r^3 cos(m theta - 5 phi)] dphi.
For the manufactured pilot m=2, epsilon=0.025, and r in [0.2, 0.8].
Fivefold periodicity reflects W7-X's five field periods; coefficients, fluxes and radial chart are *synthetic*, not W7-X measurements.

F=dA; dF=d^2 A=0. The numerical operators must implement the exterior derivative with mutually commuting coordinate differences. This tests the homogeneous Maxwell identity, not magnetohydrodynamic equilibrium, magnetic flux conservation under plasma evolution, boundary traction, the Einstein constraints, or observed W7-X performance.

## Acceptance
1. Python pilot executes on 20/40/80 angular grids with 21/41/81 radial points.
2. Max discrete dF closure error < 1e-10.
3. Manufactured analytic F error converges with order >= 1.5.
4. Record outputs and commit SHA. No change to Run 43 frozen 12/12 standards.
5. Before coupling DEE, obtain an actual VMEC equilibrium, magnetic coordinate mapping and flux normalizations. Then compare magnetic and constraint residuals on identical meshes.
6. A genuine attempt to resolve 9/12 requires exact failing tests and their discrete operators; an added potential alone cannot repair a defective inner-boundary stencil or weighted-adjoint rule.

## Reference
W7-X magnetic equilibria are commonly computed with VMEC and transformed to Boozer coordinates. Fivefold periodicity is a machine symmetry, not a substitute for actual equilibrium geometry.
