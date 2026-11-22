# Run 43 final GitHub and Physical Review D submission audit

2026-11-22 UTC · GeraldDahlbergJr/dee-toroidal-helical-omniverse

**NOT READY.** The frozen baseline is reproducible from GitHub. The later experiments support a finite-grid boundary and auditor study with explicit limitations. The exact referee report has since been recovered. The original manuscript remains missing, so direct verification of its ten numbered equations and five tables remains open. A paper restricted to the demonstrated initial-data and numerical findings can be prepared from existing evidence without new physical solves. A paper claiming unconditional numerical certification, recovery of the original gate by the experimental solvers, or general Einstein–DEE dynamics lacks essential evidence.

## 1. Repository integrity

GitHub advertises `main` at `7c39280a447f569f4791dc136bc7d25657312426`. A full clone fetched all 81 advertised branch heads at inspection time and passed `git fsck --no-dangling`. The workspace contained no pre-existing checkout to compare for unpublished changes. Detached evidence worktrees were clean after read-only replays, and the audit branch starts from the observed default-branch commit. The research refs below are pushed and remain distinct from `main`; a submission can cite pinned commits without merging experiments.

| GitHub ref | Observed tip | Role |
| --- | --- | --- |
| `archive/run43-completed-candidate-2026-10-19` | `2b4a4652e0bb57bccc3cd619f07e18ebc545f8c4` | October completion manifest and identity of the single-resolution candidate |
| `research/run43-coordinated-replay-20261103-1831` | `ff07102de88a4b63686a4f71f9fa851602ef7e0b` | Original-source replay, exact baseline files, rejected alternatives |
| `run43-final-certification-20261107` | `54183649dddf26f566b956d41618cec25b218212` | Manufactured uncertainty audit and reproducible bounds accounting |
| `research/run43-mixed-flux-sbp-sat-recovery-20261120` | `d18c986dda63c26a79642c3b606407bb339980b0` | Nodal, hybrid and conforming comparison |
| `research/run43-independent-radial-reference-20261120` | `c29162b5014de0a6ffdc33db5cc0be1304af7bfc` | Fixed-angle radial control and reference archive |
| `research/run43-angular-asymptotic-20261121` | `5996f9f5de797aba7d67fc0025b7007ec40399fb` | Angular sweep, numerical evidence introduced at `9f9d79091f5aec929ed09c56158fc3c7abcb896f` |
| `research/run43-frozen-momentum-20261121` | `a86b8f046908f852722d1d034f4e0711d5fc9466` | Term substitutions and separate coupled physical solve |
| `research/run43-inner-asymptotic-20261121` | `7368b51fb9ad4026826294837e173ec7e22f8c41` | Inner-wall Taylor and signed-deficit analysis |
| `research/run43-frozen-boundary` | `3715241456898277e290fe135222582d538c182e` | Additional frozen boundary analysis, ancestor of auditor v1 and v2 |
| `research/run43-continuum-auditor-v1` | `59823a57cf1d83c5a792f874c79e891dc6793167` | Versioned momentum evaluator and manufactured qualification |
| `research/run43-auditor-v2` | `19a0e7f94a16d759761a5fe55470263b837c9adf` | Roundoff and angular qualification introduced at `9a7dda759b7dd7f4c7d073869a1d6432a8a74eb2`; tip changes report formatting only |

The exact November 3 arrays, checkpoints, solver snapshot and reports are in the coordinated replay tree. The complete coordinated SHA-256 manifest passed from the repository root. The three manufactured certification manifests for completion, final gate and final certification also passed in the detached certification tree. The frozen diagnostic Git blob is `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22` (Git object ID). Its file SHA-256 is `36516a8d6c20699b76d814b4b5cba15378ac78948e99cf53970e33388f3b6e68`. The v2 replay checked the inherited files, the coupled arrays and the native-N reference archive and members. No evidence rewrite or branch promotion was required. [Coordinated replay and manifest](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/ff07102de88a4b63686a4f71f9fa851602ef7e0b/research/RUN43_COORDINATED_REPLAY_20261103.md)

The four Zenodo files matched their deposited MD5 checksums. The October candidate ZIP also matched the historical SHA-256 `65606d49c0f9467d40e5ad8ea4856ae572adad24f986a71abd2d7b543101bf02`. The November 2 stabilized-flux ZIP has SHA-256 `c4234f93275cabf8d00efd667b19d178212fbd0552531d263f4eee8da32ba755`. Its saved arrays independently score 12/12 with the same minimum order, but all four file hashes differ from the November 3 coordinated baseline manifest. Cite the November 3 Git commit and its repository tree `0917ebbfd9a8209081246167cce06ad64d3036e0` for the exact baseline data used here. The preparation and evidence-index ZIPs contain documentation. [Zenodo versioned record](https://doi.org/10.5281/zenodo.24330220) · [October preservation manifest](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/2b4a4652e0bb57bccc3cd619f07e18ebc545f8c4/RUN43_PRESERVATION_MANIFEST.md)

The original GitHub Actions recovery artifact, ID `14695400044`, was listed as unexpired with its recorded digest and belongs to the submitted **nodal** trial. The hybrid and 4N angular-reference arrays are absent from the inspected trees and shared workspace. Pinned generators and recorded hashes provide a regeneration route, while the public checkout supplies the hybrid summaries and native-N reference arrays. Preserve the original hybrid and 4N arrays in a durable archive with file-level hashes if those comparisons will be presented as independently replayable from archived inputs. [Nodal artifact](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/actions/runs/44872355170/artifacts/14695400044) · [v2 input status](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/19a0e7f94a16d759761a5fe55470263b837c9adf/research/run43_continuum_auditor_v2/REPORT.md)

## 2. Verified numerical results

The historical gate tests finite $p=\log_2(R_N/R_{2N})\ge1.5$ for Hamiltonian flux and order-two tensor momentum in two radial bands, over three consecutive doublings on $N=8,16,32,64$. The volume-weighted RMS uses $rq\psi^6$, with $q=R+r\cos\theta$. The inner band goes through $r=0.10$ and the outer band starts at $r=0.20$. Both wall nodes are excluded. The original data have $\psi=1$ and $W=0$ on $r=0.05,0.25$ and periodic angles. The array runner checks shape, finiteness, positivity, field hashes and exact wall values.

| Dataset and evaluator | Result | Audit evidence |
| --- | --- | --- |
| November 3 baseline, original frozen gate | **12/12 PASS**, minimum **1.7408475852068372** | Reran from the four Git-tracked NPZs. Full runner JSON matched the archived runner JSON byte for byte and all 12 checks matched the independent validation JSON. |
| Quadratic hybrid SBP–SAT, original frozen gate | **9/12 FAIL** | Recomputed each order from archived band RMS values and matched frozen acceptance JSON. Inner momentum orders 0.474713503, 0.881282789, 1.123133441. Raw-array replay was unavailable. |
| Separate coupled collocation recovery, original frozen gate | **9/12 FAIL** | Reran unchanged auditor on Git-tracked NPZs. Full JSON matched byte for byte. Inner momentum orders 0.855454968, 0.971826621, 1.127653118. |
| Hybrid arrays, compact radial second derivative plus analytic metric rule | **12/12 diagnostic** | Inner momentum 1.868352, 1.568025, 1.629676. These are results under the substituted operator, retaining the original Hamiltonian checks. |
| Other candidates | Nodal mixed-flux **6/12**, conforming quadratic control **9/12**, November 3 fourth-order **9/12**, November 3 local correction **8/12** | Archived frozen gate results for distinct configurations. |

| Baseline residual | 8→16 | 16→32 | 32→64 |
| --- | ---: | ---: | ---: |
| Inner Hamiltonian | 1.740847585 | 1.929825463 | 1.996175084 |
| Inner momentum | 1.875425682 | 1.970278594 | 1.994731881 |
| Outer Hamiltonian | 1.745663903 | 1.963989315 | 2.018376693 |
| Outer momentum | 2.229142986 | 2.132462281 | 2.071256761 |

[Baseline checks](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/ff07102de88a4b63686a4f71f9fa851602ef7e0b/research/results/coordinated-replay-20261103-1831/baseline/independent-audit.json) · [Hybrid measurements](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/19a0e7f94a16d759761a5fe55470263b837c9adf/research/run43_hybrid/measurements.json) · [Coupled report](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/a86b8f046908f852722d1d034f4e0711d5fc9466/research/run43_coupled_recovery/REPORT.md)

### Continuum auditor v2

Replayed `validate.py` at the specified commit with Python 3.12.14, NumPy 2.3.5, extended precision and one OpenBLAS thread. Its assertions passed. The JSON matched except six extended-precision-derived values on two N=64 coupled rows. The largest relative change among those six was about $1.21\times10^{-7}$. Scores, verified hashes, diagnostic classes and reported significant digits were unchanged. In particular, all 24 coupled interior outcomes are roundoff-limited and all 24 native-N reference outcomes are above the floor. Hybrid and 4N floor and screen results remain summary-only. The v2 replay compares the independent residual columns against inherited records, and does not recompute the independent quartic continuum reconstruction.

| N=64 inner band | v1 RMS | Independent reconstruction RMS | Vector gap RMS |
| --- | ---: | ---: | ---: |
| Hybrid | 1.29035e−5 | 1.15895e−5 | 7.08398e−6 |
| Coupled | 1.87209e−15 | 7.30069e−6 | 7.30069e−6 |
| Native-N reference | 8.69533e−6 | 3.15515e−6 | 9.02771e−6 |
| 4N reference | 8.98101e−6 | 1.99342e−7 | 9.02905e−6 |

For coupled fields the N=64 inner-band floor is $5.52\times10^{-13}$. The independently reconstructed last-pair inner and outer orders are 2.066 and 2.015. The shifted radial quartic stencil changes the coupled inner-band continuum residual vector by $4.44\times10^{-6}$ RMS at N=64. The shared-stencil residual measures discrete balance. The independent residual has finite-grid reconstruction uncertainty. Manufactured low-mode orders are 1.979–1.995. On the 24 manufactured cosine tests the independent angular FFT relative error was below $10^{-11}$, while the composed centered stencil is zero at Nyquist and substantially attenuates near-Nyquist curvature. Same-node FFT screening cannot expose a vanished Nyquist sine. V2's roundoff margin, angular screens and independent consistency checks are prospective diagnostics. An absolute continuum acceptance gate needs a declared tolerance and an independent residual error budget. [V2 report and measurements](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/tree/19a0e7f94a16d759761a5fe55470263b837c9adf/research/run43_continuum_auditor_v2)

### November 21 findings and supersession

- The angular sweep holds $8N$ radial intervals at each audit grid. At the finest measured angles the frozen-minus-independent vector RMS was 136–8,545 times the independent residual across the first-row and inner-band samples. At N=64 and 8N angles the inner-band frozen, independent and gap RMS values were $1.996\times10^{-4}$, $5.187\times10^{-8}$ and $1.996\times10^{-4}$. The squared operator-gap contribution was 94.98%, 92.54% and 98.83% of the three hybrid deficits under the report's signed decomposition. The smaller independent residual depends on radial refinement and join traces. These results diagnose a measured operator discrepancy. [Angular sweep](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/5996f9f5de797aba7d67fc0025b7007ec40399fb/research/run43_asymptotic/REPORT.md)
- The frozen-momentum analysis localizes the finest reference gap chiefly to the composed radial vector Laplacian (96.24% signed projection at N=64 on the 4N-angle reference). Coarse radial–poloidal metric terms substantially cancel other terms. A radial-only change can widen the coarse gap. The coupled recovery improves the finest independent inner-band momentum RMS by about 37% against hybrid, while its displacement is farther from the 4N reference and its frozen first-row RMS is $7.03\times10^{-4}$. Its frozen gate fails all three inner momentum comparisons. [Term decomposition](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/a86b8f046908f852722d1d034f4e0711d5fc9466/research/run43_frozen_momentum/REPORT.md) · [Coupled recovery](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/a86b8f046908f852722d1d034f4e0711d5fc9466/research/run43_coupled_recovery/REPORT.md)
- The inner-wall study explains the leading first-row Taylor term $hW_{rrr}/4$ in the composed radial derivative. The radial momentum contribution multiplies it by $\psi^{-10}(I+e_re_r^T/3)$. For exact smooth samples, bounded positive weights, angular steps of order $h$, a fixed-width band and a nonzero leading wall coefficient, one row of $O(h)$ error with $O(h)$ band weight gives a leading $O(h^{3/2})$ RMS. This is a truncation result. Transfer to computed fields requires control of differentiated solution error. The leading coefficient may vanish, and the finite-grid margin can have either sign. The first row contributes 98.3%, 102.8% and 97.4% of the finest signed deficit for hybrid, coupled and reference respectively. Shares above 100% reflect cancellation. [Derivation and measurements](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/tree/7368b51fb9ad4026826294837e173ec7e22f8c41/research/run43_inner_asymptotic)

The term analysis refines the earlier angular study's attribution. V2 supersedes raw order interpretations for roundoff-limited **v1 diagnostic** outcomes and qualifies v1's smooth-field consistency claim with the angular blind spot. Neither change replaces the frozen acceptance gate, baseline data or failed experiment results.

## 3. Manuscript consistency and source provenance

The follow-up recovered `Run43_Referee_Report.docx` from the versioned Zenodo review record. The Git copy at commit `9743ff740bb6a677729eb7b71cfc6b7aab335ebc` has SHA-256 `bc84f654f70b76b3d257424cf868a115070d10716fd1929388f5ee2fe48bdda7`. The verification matrix maps source evidence to all ten equations and five tables identified in the review. The original `Run43_Manuscript.docx`, corresponding PDF, LaTeX source, bibliography and supplements remain missing from the sources searched. The August preprint and older recovered DEE drafts are earlier works. [Recovery inventory](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/9743ff740bb6a677729eb7b71cfc6b7aab335ebc/research/run43_prd_reconciliation/REPORT.md) · [Verification matrix](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/9743ff740bb6a677729eb7b71cfc6b7aab335ebc/research/run43_prd_reconciliation/VERIFICATION.md) · [Review deposit](https://doi.org/10.5281/zenodo.24430854)

The archived comparison CSV confirms the referee's reconstruction correction. For the coupled N=64 inner band, the alternate stencil increases RMS from 7.30069e−6 to 7.74250e−6, or **6.05%**. The primary and alternate 32→64 orders are **2.06613** and **1.52485**. Both reconstructions share the first-row stencil, so this comparison supplies no independent first-row error estimate. [Reproduction and source checks](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/9743ff740bb6a677729eb7b71cfc6b7aab335ebc/research/run43_prd_reconciliation/REPRODUCTION.md)

| Item to reconcile against the final draft | Evidence and required treatment |
| --- | --- |
| Numbered equations and assumptions | The available replay specifies maximal conformally flat initial data, the fixed prescribed matter projections $E=(\Pi^T K\Pi+\psi^{-4}D^T KD)/2+U$, $S_i=-\Pi^T KD_i$, $G=10^{-3}$, Hamiltonian $\Delta\psi+\lvert LW\rvert^2\psi^{-7}/8+2\pi G E\psi^5=0$, and momentum $\operatorname{div}(LW)=\Delta W+\nabla(\nabla\cdot W)/3=8\pi G\psi^6S$. Match the action, conventions, chart, contractions, dimensions and normalization to a versioned source. The main-branch recovered formulation and older DEE TeX need explicit version attribution. |
| Numerical tables and convergence | Identify the dataset, operator, norm, band, grids, and evidence file for every entry. Keep finite-grid observed orders separate from analytic consistency orders, asymptotic assumptions and algebraic residuals. Use the two 9/12 scores and label any substituted 12/12 result diagnostic. |
| Boundary statements | State the two exact radial Dirichlet conditions and periodic angular conditions, the wall exclusion in band RMS, the first-interior-row defect, the strong-trace averaging convention, and the uncertainty from radial reconstruction. Dirichlet data do not prescribe the third normal derivative. |
| Scope and citations | Cite the archive DOI and exact Git commits and files above. The physical solves address the stated toroidal initial constraints in the tested configurations. Claims of unrestricted 3-D stability, dynamical recovery, full spacetime evolution or physical verification require separate evidence. Manufactured tests can exercise general 3-D inputs without extending the sector of the physical recovery solves. |
| Reproducibility and manuscript package | Include a versioned source and compiled PDF, references and equation numbering, figures/tables with traceable data, a limitations section, and a data availability statement that distinguishes DOI data from Git-only arrays. Use the recovered referee report to associate each concern with a revision or an explicitly unresolved item. Presence and adequacy in the intended draft are unverified. |

## 4. Outstanding scientific limitations

- The November manufactured certification is **NOT YET ESTABLISHED** with 0/60 total-error enclosures and 0/40 unconditional convergence comparisons. All 40 comparisons pass conditional on valid 10% endpoint error enclosures, with minimum conditional order 2.6643751474494595. The sampled Gauss-6 Jacobian inverse bounds below 1.33 apply at the stored coefficients. They do not prove an exact-assembly nonlinear neighborhood. The later interval/Taylor pilot supplied a polynomial bounding primitive and did not complete the non-polynomial integral or propagated error bounds. [Certification](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/54183649dddf26f566b956d41618cec25b218212/research/run43_final_certification/REPORT.md) · [Interval pilot](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/bbf804c9982bc20bb5dbd623d7d41142b0d5932c/research/run43_interval_taylor_pilot/REPORT.md)
- The frozen boundary discrepancy survives angular refinement and appears on the refined reference. The finite-grid evidence does not settle the limiting physical wall trace, uniform coupled coercivity or angular stability, an independently bounded absolute continuum residual, or whether an independently accurate solver can satisfy the frozen finite-grid gate. A more accurate continuum evaluator can expose the discrepancy without curing it.
- No new numerical experiment is required to substantiate the historical finite-grid gate or report the measured sensitivity as such. Unconditional certification would require rigorous fixed-field integration, arithmetic and truncation enclosures, absolute weak-defect and assembly perturbation bounds, and nonlinear error propagation to the volume and surface residuals.

## 5. Reproduction without new physical solves

From a full clone, create disposable worktrees. Use Python 3.12.14 and NumPy 2.3.5, with extended precision wider than binary64 for v2. `jq` is used for JSON comparisons. Run each shell block from the indicated directory.

```bash
git fetch origin
audit_work=$(mktemp -d)
git worktree add --detach "$audit_work/evidence" 19a0e7f94a16d759761a5fe55470263b837c9adf
git worktree add --detach "$audit_work/certification" 54183649dddf26f566b956d41618cec25b218212
export OPENBLAS_NUM_THREADS=1
cd "$audit_work/evidence"
replay=research/results/coordinated-replay-20261103-1831
sha256sum --quiet -c "$replay/SHA256SUMS"

mkdir "$audit_work/baseline" "$audit_work/coupled"
cp "$replay"/baseline/{fields_*.npz,result_*.json} "$audit_work/baseline/"
coupled=research/run43_coupled_recovery/results/fields
cp "$coupled"/{fields_*.npz,result_*.json} "$audit_work/coupled/"
python research/audit_run43_arrays.py "$audit_work/baseline"
cmp "$replay/baseline/independent-audit.json" "$audit_work/baseline/independent-audit.json"
cmp <(jq .checks "$replay/validation/baseline-independent-validation.json") \
    <(jq .checks "$audit_work/baseline/independent-audit.json")
coupled_status=0
python research/audit_run43_arrays.py "$audit_work/coupled" || coupled_status=$?
test "$coupled_status" -eq 1
cmp "$coupled/independent-audit.json" "$audit_work/coupled/independent-audit.json"

python research/run43_continuum_auditor_v2/validate.py "$audit_work/v2.json"
diff -u research/run43_continuum_auditor_v2/measurements.json "$audit_work/v2.json"
```

The coupled runner deliberately returns exit 1 for the scientific gate failure. For v2 inspect any diff in light of the extended-precision differences described above. Optional `--hybrid` and `--refined` inputs must match the recorded original hashes. The pinned reports provide full regeneration commands if an archived input is unavailable.

To rederive the hybrid score from the saved RMS values, run in the evidence worktree.

```bash
python - <<'PY'
import json
import math
from pathlib import Path

root = Path("research")
hybrid = json.loads((root / "run43_hybrid/measurements.json").read_text())
rms = {row["n"]: row["band_rms"] for row in hybrid["experiments"]["hybrid"]["levels"]}
archived = json.loads((root / "run43_frozen_momentum/results/frozen-acceptance.json").read_text())
passed = 0
for check in archived["checks"]:
    coarse, fine = check["pair"]
    band, metric = check["band"], check["metric"]
    order = math.log2(rms[coarse][band][metric] / rms[fine][band][metric])
    accept = math.isfinite(order) and order >= 1.5
    assert abs(order - check["order"]) < 1e-12 and accept == check["passed"]
    passed += accept
    print(band, metric, coarse, fine, order, accept)
assert (passed, len(archived["checks"])) == (9, 12)
PY

cd "$audit_work/certification"
for dir in run43_bspline_completion run43_final_gate run43_final_certification; do
  (cd "research/$dir" && sha256sum --quiet -c SHA256SUMS)
done
study=research/run43_final_certification
python "$study/decide.py" > "$audit_work/decision.json"
cmp "$study/decision.json" "$audit_work/decision.json"
python "$study/nonlinear/certify.py" > "$audit_work/bounds.json"
cmp "$study/nonlinear/bounds.json" "$audit_work/bounds.json"
```

## 6. PRD readiness and remaining actions

**NOT READY for the intended submission package.** The frozen numerical record can already support a limited numerical paper. The missing versioned manuscript prevents verification of the actual submission claims and of revisions addressing the recovered referee report.

1. Provide or locate the exact Run 43 manuscript source, compiled PDF, bibliography, tables and supplements. Pin their versions and complete the numbered equation, assumption, table, reference and referee cross-check against Sections 2–4.
2. Limit the paper to the verified initial-data calculations, finite-grid 12/12 baseline and diagnostic conclusions. Report both experimental 9/12 outcomes and the unresolved first-row and continuum uncertainty. Remove or substantiate any stronger claim the final text makes.
3. Publish a data availability statement that cites the exact November 3 baseline Git tree and the separate Zenodo deposit. If the manuscript relies on independently replayable hybrid and 4N comparisons, deposit the digest-matched inputs or explicitly disclose which results require regeneration and which v2 checks only use archived summaries.
4. Verify the compiled manuscript and supplements against the final citations, data links and reproduction commands before any submission.
