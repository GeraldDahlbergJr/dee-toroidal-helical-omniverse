# Run 43 operator consistency investigation

No production source, saved field, checkpoint, threshold, frozen audit or remote branch was changed. Source scripts were read from commit `12b2b2a27e09a38a304afb9dff521c814a219901` into this isolated directory; `source-hashes.json` identifies every extracted script. Existing baseline and rejected candidate conclusions remain unchanged. No publication is authorized or performed.

## Derivation from coupled constraints

The ADM constraints for maximal conformally flat initial data are R−A_ij A^ij=16πGE and D_j A^ij=8πG S^i. With γ_ij=ψ⁴δ_ij, A^ij=ψ⁻¹⁰(LW)^ij and K=0, these give

    Δψ + |LW|² ψ⁻⁷/8 + 2πGE(ψ)ψ⁵ = 0
    ∂j(LW)ij = 8πG ψ⁶ Si
    (LW)ij = ∂iWj + ∂jWi − (2/3)δij∂kWk.

For original DEE kinetic matrix KAB, fixed normal velocities Π and fixed gradients DAφ, E(ψ)=temporal/2+ψ⁻⁴ spatial/2+U and Si=−KAB ΠA DiφB. The physical residual is ψ⁻¹⁰∂j(LW)ij−8πGψ⁻⁴Si. These constraints are necessary initial-data equations; they do not solve the coupled matter evolution equations.

Commuting smooth Cartesian derivatives gives ΔW+grad(div W)/3. Discrete chart gradients Gi=Σα aiαDα do not commute and their square sum need not equal a separately discretized Laplacian L. Thus the exact identity is

    (T−E)i=(Σj GjGj−L)Wi+Σj[Gj,Gi]Wj.

The original-source fourth-order direct Hessian candidate H and frozen T additionally differ through pure derivative composition and discrete product rules. On the first interior row the frozen three-point D satisfies

    (DD−Drr)f₁=(−f₀+3f₁−3f₂+f₃)/(4h²)
               =h f'''(r₁)/4+O(h²).

Its radial vector contribution is (I+er erᵀ/3) times that defect. Dirichlet W=0 does not impose W'''=0. An O(h) residual on one row has an asymptotic fixed-band RMS contribution O(h³ᐟ²), but changing sampled derivatives and coefficients prevents a finite-pair lower bound of 1.5.

## Fresh executed evidence

`operator-comparison-fresh.json` comes from the unmodified source comparison script at n=8,16,32 and verifies the complete expanded-versus-tensor identity with relative errors 2.02e−16,2.50e−16,2.27e−16. Sampled direct 4D DEE source comparisons differ by at most 1.07e−14. Smooth tests isolate radial closure errors; these are manufactured tests only.

`verify_operator_independent.py` is a separate implementation, importing only the hash-verified frozen Audit. It forms the analytic Cartesian Hessian of W=c(.3−r−.0125/r) and continuum forcing ΔW+grad(div W)/3, then applies the unchanged frozen derivative on n=8,16,32,64. Both radial Dirichlet conditions are satisfied analytically. It independently checks the first-row matrix coefficients and tensor/commutator identity. The frozen inner M2 truncation RMS values are 5.271231524,4.252149017,2.449881290,1.139112621; orders are 0.309947903,0.795480313,1.104801455. Every order fails 1.5, even for exact continuum fields. This is an operator counterexample, not original-source physical validation or a substitute acceptance run.

## Scientific conclusion

The frozen finite-grid gate can reject exact continuum fields and does not independently rank continuum accuracy. It remains the required unchanged acceptance gate. Its defect must not be inserted into the physical forcing or used to tune solution coefficients. Direct Hessian derivatives are supported by continuum equations and accuracy testing, but existing fourth-order original-source results remain rejected 9/12. The endpoint candidate remains rejected 8/12. A fresh physical solve and all twelve unchanged frozen orders remain necessary for any new candidate.

No equation-derived new correction meeting the original-source frozen gate is established here. In particular, reducing linear tolerances cannot remove a ~1e−3 operator discrepancy when the physical algebraic solve residual is ~1e−13 (historical first-row adjudication; not presented as newly executed original-source evidence). The retained 12/12 baseline also requires supplementary independent physical validation before stronger accuracy claims.

Reproduce from this directory:

    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python source/research/compare_run43_operators.py NEW-comparison.json
    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python verify_operator_independent.py

The second command writes its specifically named evidence files; preserve this directory before replaying elsewhere. Historical sources and all original physical artifacts remain in Git and are never overwritten.

## Direct Hessian correction supported by the equations

For a scalar Cartesian component f, the orthonormal toroidal Hessian components follow directly from the flat connection (q=1+r cosθ):

    Hrr=frr
    Htt=fθθ/r²+fr/r
    Hpp=fφφ/q²+cosθ fr/q−sinθ fθ/(r q)
    Hrt=frθ/r−fθ/r²
    Hrp=frφ/q−cosθ fφ/q²
    Htp=fθφ/(r q)+sinθ fφ/q².

Rotate these tensors with the orthonormal basis er,et,ep; trace gives Δf. Then the vector equation is (HW)i=ΔWi+(1/3)Σj Hess(Wj)ij. This removes discrete derivative-composition artifacts without altering continuum sources. The existing fourth-order candidate already implements this equation-derived remedy. Its failure of the retained frozen gate cannot justify adding back the known T−H defect as a counterterm. A new investigation must supply independent accuracy/stability justification before choosing any additional discretization correction.

The fresh full sparse first-row certificate is preserved under `fresh-first-row-certificate/` (n=8,16). It reproduces relative coefficient identities 5.1641524e−16 and 5.9366461e−16, including all input boundary coefficients and all angular/Leibniz terms.

## Fresh original-source fourth-vector physical identity

Both the preserved sparse-Hessian diagnostic and a newly written independent array-only fourth-Hessian reconstruction were executed on the genuinely fresh solver output `solver/source-current/research/results/live-20261103-fourth` at n=8,16,32,64. Full row arrays are under `fresh-physical-row-identity/` and `fresh-independent-array-identity/`; field hashes are recorded by the independent implementation. These are fresh analyses of fresh solved physical fields, not manufactured solutions.

| n | Independent physical algebraic first-row RMS | Operator-mismatch first-row RMS | Frozen first-row RMS |
|---:|---:|---:|---:|
| 8 | 3.354815e−13 | 1.252558428e−3 | 1.252558427e−3 |
| 16 | 5.685446e−13 | 4.213922070e−4 | 4.213922072e−4 |
| 32 | 2.007364e−14 | 6.119317347e−4 | 6.119317347e−4 |
| 64 | 1.032487e−13 | 6.833769963e−4 | 6.833769963e−4 |

The identity M_frozen=ψ⁻¹⁰(HW−b)+ψ⁻¹⁰(TW−HW), b=8πGψ⁶S, closes across all interior rows to at most 5.15e−18 in both implementations. Independent fourth-Hessian unweighted algebraic L2 residuals are 3.93e−12,1.85e−11,2.98e−11,1.25e−11. Differences from the sparse implementation's ~1e−13 first-row residuals are small floating-order effects; both expose the same ~1e−3 mismatch. Neither analysis changes a source or acceptance condition.

This independently confirms the physical source/ψ-factor operator discrepancy, rather than merely repeating the historical adjudication. It supplies no justified additional correction: the accurate direct Hessian already follows the equations, while restoring frozen composition error would optimize for audit cancellation. The unchanged original-source twelve-check audit remains the final acceptance authority, and the fresh fourth candidate remains rejected 9/12.
