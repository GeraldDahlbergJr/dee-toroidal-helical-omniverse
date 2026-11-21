# Run 43 physical-shell mixed-flux/SBP recovery — experimental checkpoint

Date: 2026-11-20. Branch: `research/run43-mixed-flux-sbp-sat-recovery-20261120`.

## Exact scientific scope
The new `research/run43_mixed_flux_sbp_physical.py` subclasses the actual
`TensorGrid` production-experiment geometry, retains the original DEE matter
source `8*pi*G*psi^6*S`, the nonlinear Hamiltonian coupling, Fourier scalar
preconditioning and the original 8/16/32/64 physical toroidal shell. It uses
the auxiliary strain flux `sigma = L W`, with the vector momentum weak form

    integral(z dot div sigma) = - integral(grad z : sigma)

for zero-Dirichlet vector test functions and the original exterior geometry.
The discrete negative divergence is the weighted adjoint of the discrete
Cartesian gradient, with physical nodal weight `H_volume ~ r*(R+r cos theta)`
and a diagonal-norm SBP(2,1) radial derivative.

The resulting operator is **not** the same as the frozen collocated strong
tensor-divergence diagnostic. Any discrepancy is measured by the unchanged
checker. The stabilizing O(h_angle^2) fourth angular differences remain a
separately disclosed correction.

### SBP versus SAT limitation
The original solver eliminates both radial endpoint unknowns by strongly
imposing homogeneous Dirichlet values. Consequently the symmetric consistent
Dirichlet SAT boundary terms disappear on the reduced interior test/trial
space. This is a **mixed-flux SBP weak-form transplant**; it does *not* yet
implement genuine weak SAT imposition of the radial traces. A full SAT
experiment would retain boundary unknowns, add consistent, adjoint and penalty
terms in the proper weighted norm and check coercivity on the full space.
Do not label this as completed full SBP–SAT coupling.

### Unaltered acceptance
The file `reproduce/diagnose_toroidal_boundary_convergence.py` remains
immutable with Git blob SHA `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`.
The independent `research/audit_run43_arrays.py` always checks exactly
12 orders for H and M2 on both original radial bands at resolution pairs
8→16, 16→32, and 32→64. Threshold remains p≥1.5. No source fitting,
relaxed tolerances, modified physical data, or changed band definitions.

### Execute
    OPENBLAS_NUM_THREADS=1 DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier \
      DEE_EXPERIMENT=mixed-flux-sbp-physical python research/run43_mixed_flux_sbp_physical.py
    python research/audit_run43_arrays.py research/results/mixed-flux-sbp-physical

The GitHub workflow `.github/workflows/run43-mixed-flux-sbp-recovery.yml`
automates unit tests, four physical resolutions, frozen audit and evidence
upload. A failed run must be reported as failure; no promotion or conclusion
that 9/12 was fixed is authorized until the frozen output confirms it.

## Next physics step
The independent `research/run43-w7x-flux-potential-pilot-20261120` branch
holds the synthetic fivefold magnetic potential identity test. Wendelstein 7-X
microwave plasma requires an actual equilibrium, wave absorption model,
plasma response and measured data; it remains separate and has not been
coupled into this numerical recovery experiment.
