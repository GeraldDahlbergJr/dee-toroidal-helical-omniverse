# Run 43: radial endpoint accuracy and independent solution-convergence study

## Preserved baseline (DO NOT MODIFY)
- Experimental branch: `research/run43-consistent-momentum-20261102`
- Successful workflow: https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/actions/runs/41386724276
- Evidence artifact ID: `13115731417`; expires 2027-01-31; SHA-256 archive digest `7390980db3b499bd5d15438753781cb807299c008181cfb270fa01d6a907e35b`
- Frozen audit blob: `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`
- Original candidate artifact `11847550871`; expires 2027-01-11.
- Experimental reported 12/12 frozen residual gates; min order 1.7408475852. Not full physical validation.
- Preserve all source snapshots, raw arrays, checkpoints, logs, hashes, failed candidates, and metadata. No edits to baseline or production.

## Next falsifiable threshold
1. Derive the composed radial endpoint operator and its leading local truncation error analytically, including the first interior row and radial vector momentum tensor divergence.
2. Design a consistent second-order endpoint closure. State boundary assumptions, discrete commutators, and expected principal symbol. Do not tune to the twelve frozen gates.
3. Add manufactured-solution *unit tests for truncation order only*, clearly separated from original-source physical validation; test radial endpoints, near-boundary rows, angular checkerboard modes, and refinement trends.
4. Run fresh physical source-coupled 8,16,32,64 solutions under identical documented settings on a separate experiment branch. Retain full fields and checkpoint histories for failed and successful candidates.
5. Evaluate unchanged frozen H/M2 12/12 gate and independent solution-error/self-convergence in fixed physical regions. Record pointwise, RMS and boundary error norms and numerical stability; residual gate alone is insufficient.
6. Independently rehash raw arrays, solver snapshots, provenance, frozen audit blob and numerical library versions; store immutable evidence and upload artifact. Record archive retention, upload failures and a verified off-GitHub copy where technically supported.
7. If any scientific criterion fails, document the failure; do not promote, merge, alter frozen criteria, exclude rows or claim full physical validation.

## Exit criteria
- Frozen audit SHA unchanged and all 12 physical-source residual gates >= 1.5;
- demonstrable at-least-second-order endpoint truncation and solution convergence under independent checks;
- no unresolved unstable or null modes within tested scope;
- provenance and preservation verified.
This plan is not a claim that these criteria have been met.
