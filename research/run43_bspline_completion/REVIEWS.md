# Specialist review record

The integration owner ran the single qualification. The memory specialist
owned the operator changes and bounded preflights. The numerical,
mathematical, validation, and evidence reviews used separate responsibilities.

## Numerical implementation

The weak vector form uses the full Cartesian symmetric trace-free strain
and the transpose-weighted gradient. Symmetry and trace removal give the
factor 1/2 in the weak bilinear form. Removing endpoint radial basis
functions imposes the homogeneous traces of psi-1 and W. Translating a
cardinal basis across either angular seam shifts indices by n and leaves
the modulo-n coefficients invariant. Simple interior knots give C3
continuity, including at the periodic seams. Endpoint samples alone would
supply weaker evidence because the implementation reduces angles modulo
2 pi.

The vector path would also have allocated full gradient and strain tensors
after the original forcing failure. The same slab decomposition now covers
those kernels and the scalar Newton reaction term. The operator remains
matrix-free and uses the existing tensor Poisson preconditioner.

## Mathematical verification

With D_i the Cartesian derivative reconstructed from the chart reciprocal
frame, the manufactured scalar source is minus the analytic Laplacian and
nonlinear Hamiltonian load. The momentum source is minus the analytic
vector Laplacian and one third of the Cartesian gradient of divergence,
plus the original conformal source. Integration by parts with homogeneous
radial test traces gives the weak forms in the pinned protocol. The
independent strong evaluator's frame derivatives and Hessian contractions
are consistent with the chain rule. These checks leave the source equations
and Dirichlet problem unchanged.

Spline conformity and coercivity support approximation in the principal
energy norms. The required strong residual convergence, especially on the
first interior planes, is evaluated separately. CONDITIONING.md certifies
the actual exact-quadrature coupled Jacobian at the archived coefficients.

## Acceptance and evidence

Final coupled algebraic defects are recomputed at the stored coefficients.
The field errors use the analytic solution, and strong residuals are evaluated
through local power polynomials without the weak operator. The evaluator
and assembly share the symbolic forcing definition. Independent inspection
verified that definition against the analytic continuum expressions.
The residual q6/q8 comparison is an empirical quadrature diagnostic.
A verified local Jacobian inverse alone does not propagate residual and
assembly perturbations through the nonlinear problem or bound quadrature
and floating-point errors in the strong evaluator.

Historical result files and coefficient archives are inherited from the
evidence parent. New runs write to a separate directory. SHA256SUMS seals
the new evidence. The development commit updates only the handoff outside
the new evidence directory.
