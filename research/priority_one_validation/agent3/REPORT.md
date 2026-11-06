# Agent 3: equation, source and boundary audit

Checkpoint reviewed: `5ed79712e8170e1d794b358b16918094f79431f8`.
No fields, accepted archives, solver or historical evidence were changed. No DEE solve was run.

## Independently checked equations

The available model provenance is `docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md`, a versioned transcription of an author manuscript, not an independent observational justification of the action. With signature (-,+,+,+), let field kinetic matrix be C=diag(1,K_phase), normal derivative Pi=n^mu partial_mu Phi, and spatial derivative covector D_i Phi. Varying the action gives T_munu=C_AB partial_mu Phi_A partial_nu Phi_B+g_munu L. Projection gives

- E=(Pi^T C Pi + gamma^ij D_i^T C D_j)/2+V;
- S_i=-Pi^T C D_i;
- S_ij=D_i^T C D_j+gamma_ij[(Pi^T C Pi-gamma^kl D_k^T C D_l)/2-V].

`reproduce/dee_stress_energy_3p1.py:projections` implements these expressions, including the symmetric mixed phase coupling. `check_equations.py` independently constructs a 4D metric with random nonzero shift and lapse, forms T directly, and contracts with the unit normal. Across 512 deterministic SPD-metric samples, maximum absolute discrepancies (E,S_i,S_ij,trace) are (7.37e-14,3.29e-14,1.24e-12,2.93e-13), consistent with floating-point cancellation. This extends the old zero-shift check; it tests algebra, not PDE accuracy. No threshold was fitted or changed.

For gamma_ij=psi^4 delta_ij, K=0, A_phys^ij=psi^-10 a^ij, a=LW, conformal geometry gives R3=-8 psi^-5 Delta psi and A_phys,ij A_phys^ij=psi^-12 |a|^2. Hence the physical Hamiltonian H=R3-|A_phys|^2-16piG E is zero iff

    Delta psi + psi^-7 |LW|^2/8 + 2piG E psi^5 = 0.

Tracefree conformal tensor divergence satisfies D_j A_phys^ij=psi^-10 partial_j a^ij. Because S_phys^i=psi^-4 delta^ij S_j, momentum M^i=0 is equivalent to

    partial_j (LW)^ij = 8piG psi^6 delta^ij S_j.

These are the production target equations in `reproduce/solve_toroidal_helical_coupled_constraints.py:energy,longitudinal,solve_vector,solve,diagnostics`. Tensor contraction sums all nine Cartesian entries, so off-diagonal terms count twice correctly. Cartesian components differentiate as Cartesian scalars on the chart; inserting toroidal component connection terms on top would double-count coordinate rotation. The expanded identity div LW=Delta W+(1/3)grad div W is a continuum identity; it is not an identity between arbitrary discrete derivatives. Thus a finite-grid mismatch is not evidence of a missing continuum term.

`nonlinear_initial_data.py:matter_momentum_density` uses Pi=-n·partial Phi and the corresponding positive source sign explicitly. That is a different declared convention, not a confirmed production sign error. It must not be substituted into Run43 with the same numerical Pi values without conversion.

## Quantitative normalization and geometry

Using the source seed rho^2=.945, P=(.12,-.04), d=(.25,.20), direct coefficient calculation gives

    K_phase=[[1.2835,.0945],[.0945,1.189]],
    eig(K_phase)=(1.1305957880631348,1.3419042119368652),
    E(psi)=.0127638+.068614375 (r^-2+9q^-2) psi^-4,
    S_cart=-.030316 (e_theta/r-3 e_phi/q), q=1+r cos(theta).

These agree with the production source, including G=.001 and the covector-to-Cartesian reciprocal-basis conversion. Holding E fixed when changing psi would alter the specified model; holding Pi and coordinate covectors fixed is the stated physical prescription. No confirmed normalization, sign, conformal exponent or mixed-sector error was found.

The chart embedding X=((1+r cos theta)cos phi,(1+r cos theta)sin phi,r sin theta) yields metric diag(1,r^2,q^2), volume rq dr dtheta dphi, gradient e_r partial_r+e_theta partial_theta/r+e_phi partial_phi/q. On r in [.05,.25], q>=.75 and r>0, so there is no axis/chart degeneracy. The spatial metric is positive for psi>0; this alone is not a full spacetime evolution test.

## Boundary interpretation and direct archived checks

The prescribed problem is a finite annulus with psi=1,W_cartesian=0 at both radial surfaces, periodic angular traces. Outward normals are -e_r at the inner surface and +e_r at the outer. These are mathematical Dirichlet benchmark data. They do not specify exterior matching, asymptotic flatness, traction, axis regularity or a vanishing third normal derivative. There is no action-derived authority to add an endpoint derivative condition to force a discrete audit to agree.

Fresh read-only checks of all four coordinated-replay baseline arrays (n=8,16,32,64) find finite fields, positive psi, exactly zero nodal Dirichlet errors for psi and all W components. The finest psi range is [1,1.0001670221573893]. File SHA256 values and all numbers are in `equation_checks.json`. Nodal Dirichlet agreement does not establish accurate one-sided normal derivatives, continuum endpoint constraint satisfaction, or angular seam derivative convergence. Angles store unique periodic nodes, so a duplicated-endpoint equality check is not meaningful; derivative/reconstruction periodicity needs separate evidence.

## Scope and verdict

**PASS, scoped:** action-to-projection algebra, source normalization and conformal continuum equations are consistent with the implemented model; stored nodal Dirichlet data are exact.

**NOT YET ESTABLISHED:** independent continuum H and all three M components with accepted uncertainty and p>=1.5 at the inner boundary. Algebraic correctness does not remedy a discretization/evaluator failure. This investigation identifies no new specific equation/source/boundary correction that justifies a DEE solve.

The complete initial-constraint system can legitimately be tested without solving the matter evolution equations; it cannot then establish the full Einstein–matter evolution, matter-current conservation over time, exterior physics, general nonaxisymmetric behavior or global phase winding admissibility. The seed is axisymmetric in stress although Cartesian components vary with phi. Globally admissible phase values/winding require a declared target-space convention; gradients alone suffice for this local constraint benchmark. These are limits on physical claims, not demonstrated causes of the numerical residual discrepancy.

## Reproduction

From the repository root, run `python research/priority_one_validation/agent3/check_equations.py`. It reads the immutable baseline arrays, performs 512 action-projection samples, and writes only `research/priority_one_validation/agent3/equation_checks.json`. It does not call a PDE solver. NumPy is sufficient for this script; Python's import of the source projection module performs no computation or archive writes.
