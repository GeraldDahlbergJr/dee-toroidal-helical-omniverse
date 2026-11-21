# Run 43 independent radial reference

The reference residual accounts for 95.33%, 91.16% and 97.77% of the hybrid's squared-RMS shortfall on the three inner momentum gates. The shares include the sign of each contribution. The verified hybrid score remains 9/12. The native quartic discretization scores 8/12 at the original resolutions.

## Reference problem

Base commit `d18c986dda63c26a79642c3b606407bb339980b0`. The reference uses continuous quartic Lagrange elements on uniform radial nodes and six-point Gauss integration. It shares no radial differentiation, penalty, interface unknown or scalar radial flux routine with the hybrid solver. At a given angular N it retains the hybrid's centered theta derivatives, equivariant azimuthal derivative, directional damping coefficients and scalar theta face flux. Each radial sequence uses N, 2N, 4N and 8N radial intervals. The original audit nodes coincide with quartic nodal values at every refinement.

With J=r(R+r cos θ), W is obtained from the positive bilinear form (Lv,LW) weighted by J/2, plus the hybrid angular stabilization, equated to the negative original momentum body load. The boundary work before trace restriction is the pairing of v with (LW)n in surface measure J dθ dφ. Continuity cancels interface traction terms in the weak form, and the physical Dirichlet test trace vanishes. The vector stiffness is a sum of strain and damping squares. Within this equivariant sector, zero damping removes horizontal displacement and theta dependence, zero strain forces constant vertical displacement, and the wall conditions eliminate the constant. This gives finite-grid positivity. No uniform coupled stability estimate was obtained.

The scalar solve uses the centered angular face-flux operator at quadrature points and independent weak radial differentiation. Its wall values are one. The vector source scales as ψ⁶ and the Hamiltonian includes the locally evaluated strain squared, kinetic matter, and potential. A Newton solve for the scalar field alternates with the linear vector solve. All 16 runs converged in two outer iterations. The largest final assembled vector L2 residual was 1.55e−16 and the largest scalar residual in the weak equations was 1.92e−17. These depend on quadrature normalization.

On the native grids, comparison to the original DEE source at 40 quadrature points per grid gave maximum absolute disagreement 7.11e−15. The radial discrete Green identity EᵀQD+DᵀQE=B had maximum defect 1.12e−15. Here E evaluates the continuous nodal trial space at quadrature points, D evaluates its radial derivative and Q contains quadrature weights. The vector stiffness symmetry defect was at most 2.14e−14.

## Verification and uncertainty

The manufactured displacement is vertical, with amplitude exp(x)x(1−x)+sin(2πx)/2+x(1−x), x=(r−a)/(b−a). Both wall traces vanish and its inner radial curvature is −50. For this radial verification the force uses analytic radial flux differentiation and the exact prescribed semidiscrete angular adjoint, implemented separately from the assembled stiffness. It therefore measures radial accuracy at fixed angular N. For N=8 the maximum nodal displacement error decreased from 4.66e−4 to 8.11e−8 across the four radial refinements. For N=64 it decreased from 8.11e−8 to 4.83e−12.

The uncertainty indicator is the hybrid-weighted RMS of the pointwise change in the reference residual on the original nodes between 4N and 8N radial intervals. The preceding changes, field changes and manufactured errors are in refinement.csv. The indicator is an empirical radial sensitivity; it does not certify the reference error or constrain angular error. The two finest angular grids approach floating-point sensitivity. A termwise perturbation of a reference energy E=R² has magnitude at most 2Rε+ε² for a perturbation of measured norm ε. Applying this to the last refinement on both members of each pair gives reference gate-deficit sensitivities of 4.14e−13, 3.86e−15 and 3.81e−17, respectively. The corresponding reference deficits are 4.80e−7, 1.03e−7 and 1.64e−8.

## Residual accounting

Each archived hybrid field was verified against its hash at the base commit. On each audit grid let H and F denote the frozen momentum residual vectors evaluated on the hybrid and finest radial reference fields. Define D=H−F point by point. For every band and row the same normalized weights J ψ_hybrid⁶ give

    ||H||² = ||F||² + ||D||² + 2<F,D>.

The rows.csv file includes every one of the 30 inner-band rows, its radius and normalized band weight, its row-local decomposition, the last radial change, and its additive contributions to the band energies. The largest band closure error is 1.39e−22. A negative cross term means partial cancellation. A residual ratio alone cannot distinguish cancellation from a smaller underlying error.

For a transition, define the signed gate deficit for each squared-energy term as its fine value minus one eighth of its coarse value. A positive total means the p≥1.5 gate fails. Linearity gives a signed contribution from each term. The difference-only shares of the three deficits were +0.11%, −0.28%, −0.17%. Cross-term shares were +4.56%, +9.11%, +2.41%. Thus the reference term controls every pair, including the coarsest transition where the difference-only order is also below threshold.

These results isolate the effect of changing the radial solver while holding the angular formulation fixed. The reference term includes frozen differentiation error, finite-angular-resolution error in its underlying fields, and their interaction. Separating those effects requires an independent angular refinement control. The data support retaining the hybrid boundary formulation for now.

## Evidence and reproduction

- [Full band decomposition](results/bands.csv), [first row and every other inner-band row](results/rows.csv), and [signed pair deficits](results/pairs.csv).
- [Radial refinement and manufactured errors](results/refinement.csv), [16 solver histories and array hashes](results/histories.json), [source and dependency provenance](results/provenance.json), and [source and Green-identity checks](results/verification.json).
- [Unchanged native-grid auditor output](results/native-audit.json). The four failures were the first inner Hamiltonian pair and all three inner momentum pairs.
- [Reference radial fields](results/Run43_Radial_Reference_Fields.zip). SHA-256 `8dcd1904f098b1e0748078f8fea63ce009c253dcb7f06b75be70004f3f97b4fa`. Every reduced field archive contains ψ and the three cylindrical-frame displacement components on the full radial solve grid. The comparison code rotates them back to Cartesian components on the original audit nodes.

From the repository root with Python 3.12, NumPy 2.3.5 and SciPy 1.17.0, use fresh directories. The hybrid directory can also be taken directly from the original Run43_Recovery_Evidence.zip. compare.py checks its field hashes before computing the decomposition. The native audit runs on N radial intervals for each N and writes all 12 checks; finer references are used only for attribution.

~~~bash
export OPENBLAS_NUM_THREADS=1
python research/run43_hybrid/physical.py hybrid
python research/run43_reference/reference.py radial --factors 1 2 4 8
python research/run43_reference/compare.py radial hybrid comparison
python research/run43_reference/verify.py comparison/verification.json
~~~
