# Nonlinear toroidal initial-data benchmark

This extends the immutable Run 38 source baseline (2ff3b7e) using a conformally flat spatial metric, maximal slicing K=0, and a longitudinal tracefree extrinsic curvature from a three-component Cartesian vector potential W. It solves the Hamiltonian and all three momentum equations together on the shell r=0.05 to 0.25, R=1. Angular boundaries are periodic; radial benchmark boundaries are psi=1 and W=0. G=0.001 is an explicit numerical coupling, not an inferred physical unit calibration.

With Cartesian background components, the equations are

- Delta psi + (LW)^2 psi^-7 / 8 + 2 pi G E(psi) psi^5 = 0.
- Delta W + grad(div W)/3 = 8 pi G psi^6 S_background.

The physical field normal derivatives and coordinate gradients remain fixed at the Run 38 values. The energy is recalculated with the physical inverse metric: E(psi)=E_temporal+U+psi^-4 E_gradient. The mixed phase kinetic term is retained. The frozen source stress is phi-independent, even though the underlying phase is theta-3 phi. Consequently this is not a validation of generic phi-dependent matter.

The old Hamiltonian precursor used nested centered radial first derivatives with zero endpoint derivatives, and held energy fixed under metric changes. The corrected version uses the expanded toroidal Laplacian with nearest-neighbor second differences, includes actual Dirichlet boundaries, and recomputes matter energy.

## Reproduction

Install numpy and scipy, then from the repository root run:

```
OPENBLAS_NUM_THREADS=1 python reproduce/solve_toroidal_helical_full_constraints.py
python reproduce/solve_toroidal_helical_nonlinear_constraints.py
python -m unittest discover -s tests -q
```

The coupled script evaluates grids (9,16,16), (17,32,32), (33,64,64). An FFT/sine inverse Laplacian preconditions Newton-Krylov. It independently evaluates a conservative scalar face-flux Hamiltonian and differentiates the constructed LW tensor for momentum, instead of simply reporting the nonlinear solver residual. It also checks the scalar and vector operators using manufactured Cartesian quadratic fields that vary in all three chart coordinates. Momentum errors are reported both across the shell and on the fixed physical band r=0.10 to 0.20; the two-cell diagnostic alone does not use a fixed physical band.

## Result and limits

All three nonlinear solves converge with positive finite psi. The finest solve yields physical Hamiltonian max error about 1.25e-10 and conformal momentum equation errors below 2.4e-16. Independent full-shell Hamiltonian RMS remains about 6.15e-5; independent physical momentum RMS in the matched band is approximately (9.60e-6,9.60e-6,1.25e-5). These decrease under refinement, but cannot be replaced by the far smaller solver residuals. Manufactured scalar/vector errors approach second-order refinement. Boundary behavior remains to be assessed more fully.

The artifact status is SOLVER_PASS_VALIDATION_PENDING. This does not certify a generic 3-D constraint-satisfying continuum solution, physical boundary conditions, coupled evolution, or experimental validity. The next validation should assess boundary convergence and a genuinely phi-dependent stress profile while keeping the frozen baseline intact.
