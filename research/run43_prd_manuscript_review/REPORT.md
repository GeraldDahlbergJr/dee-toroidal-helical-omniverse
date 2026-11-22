# Run 43 PRD manuscript-to-evidence reconciliation (22 Nov 2026)

**Status: REVIEW DRAFT, NOT YET READY FOR FINAL APS SUBMISSION.** The specific 11-page manuscript `Run43_PRD_Submission.pdf` was inspected against frozen research reports, relevant source code, and archived measurements. This review is documentary, not a new physical-solver replay.

## Ten numbered equations

- (1)-(2): kinetic metric, potential, and 3+1 projections align with `reproduce/dee_stress_energy_3p1.py`.
- (3)-(5): Hamiltonian, vector momentum and weighted RMS acceptance definition align with the prescribed-source toroidal solver and frozen auditor.
- (6): exact rational Taylor moments of `(3f0-5f1+f2+f3)/(4h^2)` at row 1 yield `f'' + h f'''/4 + 5 h^2 f''''/24 + O(h^3)`.
- (7): the fixed-smooth-field weighted-band squared error is `C^2 h^3 + O(h^4)`; computed-field transfer needs differentiated-error control. This is not a solver-convergence theorem.
- (8): the v1 Hessian/metric divergence matches `research/run43_continuum_auditor_v1/auditor.py`; the original PDF contained broken math delimiters requiring typesetting correction.
- (9): v2's empirical uncancelled roundoff scale agrees with its archived report, but has no universal verified guarantee.
- (10): centered-first-difference composition has response `[sin(kh)/(kh)]^2` and an exact Nyquist cosine null.

## Five numbered tables

1. Frozen gate: baseline **12/12**, hybrid **9/12**, coupled **9/12**. All three failures in each experimental variant are inner-band momentum pairs.
2. First-row signed allocations in `research/run43_frozen_boundary/results/deficits.csv`: hybrid **7.0/69.7/108.3%**, coupled **10.6/81.8/112.8%**, with cancellation.
3. v1 versus independent reconstructed field RMS agrees to displayed precision with v1 and coupled comparison records; coupled v1 inner RMS at N=64 is floor-limited, not a continuum error bound.
4. N=64 vector-gap records agree with v2: hybrid **7.08398e-6**, coupled **7.30069e-6**, N reference **9.02771e-6**, 4N reference **9.02905e-6**.
5. v1 manufactured last-pair orders and v2 angular-mode response match archived reports. These do not establish full coupled stability.

## Source narrative inconsistency resolved

The primary N=64 coupled independent inner residual is **7.300687932e-6**; the alternate local quartic is **7.742500202e-6** (a **6.05% increase**), and their vector RMS delta is **4.443972621e-6**. These numbers are verified in `research/run43_coupled_recovery/results/comparison.csv`. The same report's prose saying a 33% *decrease* at N=64 is inconsistent with its CSV. Use CSV measurements and treat the comparison only as reconstruction sensitivity.

## Submission corrections / remaining gates

- A corrected local 11-page review-draft PDF fixes the broken display of (8) and reconciles the manuscript AI statement with the saved APS disclosure, naming ChatGPT GPT-6 and Codex. It has **not** been uploaded to APS.
- APS draft `es2026nov22_675` remains **Awaiting submission**. Do not press final Submit without approval.
- Data-availability wording must distinguish public code and summaries from missing hybrid/4N raw arrays. The independent raw-array replay and rigorous error bounds are not complete.
- The pinned audit `research/run43_prd_submission_audit/REPORT.md` could not see the exact manuscript. This file closes the **documentary comparison** only, not a new numerical replay or formal referee resolution.
- No claim is made of general 3-D nonlinear stability, Einstein–matter evolution, or empirical physics prediction.

Relevant pinned refs: baseline/coordinated `ff07102de88a4b63686a4f71f9fa851602ef7e0b`; frozen boundary `3715241456898277e290fe135222582d538c182e`; coupled `a86b8f046908f852722d1d034f4e0711d5fc9466`; v1 `59823a57cf1d83c5a792f874c79e891dc6793167`; v2 `19a0e7f94a16d759761a5fe55470263b837c9adf`; prior PRD audit `9e9d69b9a6027d685ceba406a728b785d92b90ee`.
