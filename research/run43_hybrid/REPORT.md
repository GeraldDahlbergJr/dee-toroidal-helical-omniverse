# Run 43 physical-source recovery

The submitted nodal SBP(2,1) mixed-flux solver scored **6/12** on the original four grids. A quadratic hybridized SBP–SAT trial scored **9/12**. The conforming control on identical quadratic elements also scored **9/12**. Each score comes from a separate invocation of the unchanged array auditor. The validated baseline remains a separate result.

## Provenance

[Original workflow](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/actions/runs/44872355170) at commit `1086093db093359f4154658f1536fb110b22f5ef` completed the four solves, uploaded the arrays and failed its final scientific gate. The audit step used continue-on-error. Its recorded outcome was failure. The [original archive](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/actions/runs/44872355170/artifacts/14695400044) has SHA-256 `e10904e83e086dab1cac372361066918b32c3779adb069ca4108b6515dddd42e`. Local re-evaluation verified every field-file hash and boundary trace. The frozen diagnostic Git blob is `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`, and its array-audit wrapper is `a9505aeb2f958e445728e1c5336e464ea6f93963`. Full-precision orders, field and code hashes, residuals and diagnostics accompany this report in measurements.json.

## Discrete operators and wall terms

For the tested nodal transplant, let M be the diagonal radial SBP norm times the angular quadrature and the toroidal Jacobian J = r q, with q = R + r cos θ. Write L for the discrete Cartesian trace-free strain and use the Frobenius inner product. The negative momentum stiffness before damping has bilinear form −(Lv,Lw)_M/2 on the space of zero wall traces. The boundary work before trace elimination is the surface pairing of v with (Lw)n, using n = ±e_r and surface weight r q dθ dφ. Algebraically eliminating wall values suppresses this term. The pure flux satisfied the weighted adjoint pairing to 1.6e−14 relative or better in the recorded probes. The original θ damping instead multiplies a fourth difference by a variable weight outside the differences. It breaks weighted symmetry of the complete operator.

Even in one dimension with constant coefficients, the first interior row of −H⁻¹DᵀHD for SBP(2,1) is (u₃−3u₁+2u₀)/(4h²). It gives 3/4 on u=x²/2. An order-one error on one row contributes order h^0.5 to a fixed-band RMS with bounded positive weights. The full-shell manufactured test with W = k²(r−a)(b−r)e_z, k = π/(b−a), exposes the defect with nonzero wall curvature. From N=8 to 64 the middle-row RMS falls from 1.628 to 0.02543 while the first-row RMS changes from 121.0 to 137.1. At 64 the actual physical first row contributes 99.91% of the squared weighted inner residual.

The second experiment uses quadratic elements spanning pairs of the original radial cells. The three original element nodes have Simpson weights h(1,4,1)/3 and the quadratic Lagrange derivative D. They obey HD+DᵀH=diag(−1,0,1). Both element endpoints carry independent volume values. A separate face displacement λ provides a single-valued physical trace. The physical wall data prescribe λ=0. The volume displacement at a wall remains an unknown. For j=u−λ and t=(Lu)n, the positive stiffness is

    a(u,v) = (Lu,Lv)_M/2 − <t(u),j(v)> − <j(u),t(v)> + <τ j(u),j(v)> + damping.

Face brackets carry surface weight J. The negative physical body force loads the volume equations. The face equations enforce balance of numerical traction. The penalty τ=8/w, w=h/3, follows an energy bound independent of the acceptance data. At each endpoint, 2|t·j| ≤ (w/4)|Lu|²+(4/w)|j|². Thus the quadratic form retains positive strain and jump terms. Angular damping is a weighted square of the undivided periodic second difference scaled by angular spacing, with coefficient 0.25 and the original directional metric factors. The conforming control restricts this exact stiffness to j=0. This isolates the effect of the SAT coupling. Shared quadrature and angular damping also differ from the nodal transplant.

The physical source is invariant under azimuthal rotations when the Cartesian vectors rotate with the frame. The new solver uses that invariant sector. Scalars are independent of φ, and the Cartesian vector rotates with φ. Centered differentiation of that vector gives sin(Δφ)/Δφ times the rotation generator. The full-grid strain inner product on this sector reduces to the meridional weighted product times 2π; the omitted common factor cancels. The scalar solver reuses the original face-flux stencil and physical matter projections. It averages squared element strain at shared scalar nodes. Audited Cartesian arrays contain the simultaneous physical solution's primary face trace at element interfaces and its element midpoint values elsewhere. Raw broken displacement and trace unknowns are archived separately. No field smoothing, boundary clipping or audit-specific extrapolation is used. These trials validate only the equivariant sector.

## Diagnostics and limits

The zero-energy conditions on the represented sector eliminate the displacement kernel. Vanishing angular damping enforces θ-independent displacement and zero horizontal components, vanishing strain makes the vertical displacement constant on each element, and vanishing jumps with zero physical wall traces forces those constants to zero. This establishes finite-grid positivity of the vector stiffness on that sector. A uniform coupled coercivity or nonlinear stability estimate remains open.

The hybrid stiffness symmetry defect was at most 7.11e−15 in maximum entry norm. The diagonally scaled eigenvalues nearest zero were positive at 0.002930, 0.0007583, 0.0001912 and 0.00004791. The manufactured solution W=sin(k(r−a))e_z uses the analytic force ΔW+grad(div W)/3. Its nodal displacement RMS fell from 0.006798 to 0.0001163. The same test separates the frozen differentiation defect on analytic samples from the twice-differentiated numerical solution error. The respective inner RMS values at 64 were 0.9481 and 0.4101, with total 0.9988. Their norms do not add linearly. This decomposition concerns the manufactured field. The exact physical-source solution is unavailable.

At N=64, the hybrid physical algebraic vector residual was 9.91e−17, the scalar interior residual maximum was 5.41e−11, and the largest displacement jump was 6.79e−10. Physical inner momentum RMS fell from 0.00329945 for the nodal transplant to 0.000203148 for the hybrid. Its first interior row still contributed 51.31% of the squared inner-band residual. The matched strong control recovered the same three outer checks. These data do not identify a unique leading coefficient in the inner physical residual or establish that weak wall enforcement caused the recovery. The original transplant's fixed-band behavior has a measured first-row mechanism. The remaining inner failure of the quadratic method requires separating physical-field discretization error from the frozen differentiation defect.

## Reproduction

Run from repository root with Python 3.12, NumPy 2.3.5, SciPy 1.17.0 and mpmath 1.3.0. The original workflow used Python 3.12.16. The follow-up used 3.12.14. Both trial drivers require fresh output directories.

~~~bash
export OPENBLAS_NUM_THREADS=1
python research/run43_hybrid/physical.py output-hybrid
python research/run43_hybrid/physical.py --conforming output-control
python research/run43_hybrid/row_diagnostics.py output-hybrid
python research/run43_hybrid/strong_diagnostics.py strong-probe.json
python research/audit_run43_arrays.py output-hybrid
python research/audit_run43_arrays.py output-control
~~~

Each audit writes independent-audit.json and returns exit code 1 for the observed scientific failures. The independent auditor is a separate computation using the original diagnostic implementation. External scientific replication remains open.

## Next experiment

Use an independently discretized smooth radial reference for the same physical forcing and angular operators. Evaluate it at the original four grids, compare the inner row profiles, and separate the frozen differentiation defect from the displacement error after differentiation. Use a radial-only refinement control for attribution outside the twelve acceptance checks. Retain the original four paired grids for all formal scores. Compare traction residuals on broken and conforming elements before proposing another boundary change. The electromagnetic pilot remains separate.
