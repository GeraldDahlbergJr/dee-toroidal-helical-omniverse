# Einstein–DEE new physical test — specification decision record v0.1

Parent evidence: `21683875ed93cd413d3b4472ba360c6e21a70a9d`.
Parent correction/protocol: `0e635ed9fea802398c23d2826d51574f3425dfd1`.
Status: **DESIGN / NOT EXECUTABLE**. Scientific verdict: **NOT YET ESTABLISHED**.
This is a NEW experiment, not reconstruction of the historical Gaussian run. No PDE solves or physical tests are authorized.

## Decisions adopted for the proposed test design
- Target: a nontrivial, stationary, axisymmetric, z-translation-invariant weak-field Einstein–DEE matter configuration in a finite cylindrical annulus, as proposed in `NEW_PHYSICAL_TEST_PROTOCOL.md`. This tests local annular consistency only, not a complete toroidal wormhole, global spacetime or nonlinear stability.
- Coordinate convention: `(t,r,phi,z)`; signature `(-,+,+,+)`; `c=1`; flat cylindrical background `diag(-1,1,r^2,1)`. Evaluate all ten Einstein components and all four harmonic-gauge components, with background-frame normalization.
- Matter: original three-scalar action and its Euler–Lagrange equations, not an arbitrary Gaussian-multiplied stress tensor. Use the polynomial potential/kinetic realization only if its coefficients and physical interpretation are explicitly fixed. Nonzero winding must respect phase topology.
- Metric: retain ten components until justified symmetry/gauge reduction; never infer a complete Einstein solution from the isolated phi-z twist equation.
- Use the existing protocol's preregistered error thresholds provisionally: first-order normalized component residual <=1e-6, weak-field frame norm <=1e-3, and separately bounded numerical uncertainty <=10% of each acceptance threshold. These are proposed methodological criteria, not measured results.
- Validation sequence after separate execution approval: analytic/implementation controls; matter closure and conserved source; gauge and all component Einstein residuals; compatible boundary matching; independent verification. Controls never count as a DEE physical PASS.
- Keep original manuscripts, historical PASS, and all frozen Run 43 paths unchanged. The Eq. (34) formal correction stands.

## Explicitly unresolved: must freeze before execution
1. **Physical parameter table and units**: length L, G in consistent units, rho and phase scales, lambda_V, v, a, b, c_mix, omega_Theta/Psi, m_Theta/Psi, k_Theta/Psi; select a physically motivated nontrivial family, not convenient post-hoc values.
2. **Weak-field family**: choose between scaling the full matter action by epsilon (a new gravitational-coupling assumption) and a derived field/coupling family. Neither has yet been approved. Specify metric/matter order and epsilon sequence.
3. **Radial domain and boundary compatibility**: r_in, R, interior/exterior continuation or explicit supporting stresses, gauge-compatible gravitational data, matter radial values/currents, junction conditions, and proof the radial constrained system is not overdetermined.
4. **Numerical implementation**: algorithm, base N, precision, quadrature/differentiation, tolerances, fixed dimensional normalization scales, certified error bounds, and exact run manifest.
5. **Independent verification**: verifier implementation, reviewer availability, version/commit locks and evidence output paths.

## Go/no-go
- **GO**: read-only derivation, analytical boundary-condition counting, parameter dimensional analysis, source-family feasibility design, and documentation of the above choices.
- **NO-GO**: PDE solves, physical residual evaluations, numerical certification, or declaration of Einstein–DEE physical consistency.
- Do not select arbitrary physical constants or fabricate a boundary continuation merely to make the manifest appear complete. A feasible nontrivial source and compatible boundary data are prerequisites.
- This record is an actionable specification checklist, **not** a claim that the test is fully specified. Obtain explicit author approval for a frozen executable manifest before running the experiment.
