# Run 43: composed endpoint derivative error budget

The saved fourth-order vector candidate's three inner-M2 failures are explained
quantitatively by a slowly converging endpoint contribution plus grid-dependent
cancellation with the rest of the frozen residual. The leading Taylor term alone
is not quantitatively sufficient at these resolutions. No solver improvement
satisfying both the frozen gate and supplementary accuracy is established here;
no new candidate is promoted.

Analysis branch: `research/run43-endpoint-error-budget-20261103`, based on
`188cbd41bcfb559c62f030c9d329aa11b8a1ac5c`. The prior restoration/recheck is used,
without rerunning restoration, baseline verification, or nonlinear solves.
Only saved fourth-vector n=8,16,32,64 arrays are evaluated. The original November
2 baseline remains accepted at 12/12; the fourth-vector candidate remains 9/12.

## Derivation and exact isolation

Let x=r-rmin and h=.2/n. The frozen three-point endpoint derivative has
D f(0)=f'(0)-h^2 f'''(0)/3+O(h^3), while the centered interior derivative has
D f(x)=f'(x)+h^2 f'''(x)/6+O(h^4). At the first interior row,

    (D D f)_1 = (3 f_0 - 5 f_1 + f_2 + f_3)/(4 h^2)
               = f''(r_1) + h f'''(r_1)/4
                 + 5 h^2 f''''(r_1)/24 + O(h^3).

The sign is positive. The mismatch of endpoint versus centered derivative error
is -h^2 f'''/2, multiplied by the outer centered coefficient -1/(2h).

For isolation only, use the existing matched endpoint derivative with weights
[-2, 7/2, -2, 1/2]/h. Its leading error equals the centered error. Keep all
interior and angular stencils fixed. Let F be the original frozen physical M2
residual and B the supplementary matched residual on the SAME saved field.
The exact inner-band difference E=F-B is supported only at row 1:

    E_1 = psi_1^(-10) (h/4) (I + er er^T/3) delta3 W,
    delta3 W = (-W_0 + 3W_1 - 3W_2 + W_3)/h^3.

The radial basis er is independent of r. The projection comes from the radial
principal symbol of Delta W + grad(div W)/3. This identity includes the physical
psi factor; the unchanged matter source cancels. The outer closure difference
is outside the inner band. The exact identity was checked on every saved field.
No modified diagnostic is used for acceptance.

For comparison with the Taylor term, estimate W''' at row 1 using polynomial
stencils of 5,7,8 points. The 8-point estimate gives relative vector differences
from exact E of 42.21%, 30.29%, 18.87%, 10.60%, respectively; directional cosines
exceed .99995. Seven versus eight points agree much more closely than either
agrees with E: the difference largely reflects finite-h higher terms, not merely
third-derivative estimator sensitivity. This is an estimate on discrete fields,
not an independent exact continuum derivative. The exact delta3 identity avoids
that estimation uncertainty entirely.

## Measured inner-band budget

All norms and inner products below use the frozen band r in (.05,.10], including
its .10 row, and the frozen volume weights r*q*psi^6. There are 2,4,8,16 radial
samples. E and B are vectors, so their norms must not be added as scalars.

| n | Frozen ||F|| | Endpoint ||E|| | Remaining ||B|| | cos(E,B) | Cancellation fraction |
|---:|---:|---:|---:|---:|---:|
| 8 | 9.839653e-4 | 1.020051e-3 | 1.880198e-3 | -.940473 | .788404 |
| 16 | 7.824655e-4 | 7.578015e-4 | 1.150057e-3 | -.737017 | .677234 |
| 32 | 4.324281e-4 | 4.311772e-4 | 4.684305e-4 | -.540525 | .538675 |
| 64 | 1.999734e-4 | 2.009109e-4 | 1.527242e-4 | -.386204 | .372124 |

Cancellation fraction is -2<E,B>/(||E||^2+||B||^2). The exact identity
||F||^2=||E||^2+||B||^2+2<E,B> closes to floating-point accuracy. At n=64 its
three terms are 4.0365191e-8 + 2.3324693e-8 - 2.3700526e-8 = 3.9989358e-8.
The endpoint term alone is slightly larger than F; that does not mean it accounts
for more than 100% of an additive error budget. Cross terms are essential.

| Pair | Endpoint order | Remaining order | Order ignoring cross term | Cancellation order shift | Actual frozen order |
|---|---:|---:|---:|---:|---:|
| 8→16 | .428749 | .709179 | .635168 | -.304588 | .330580 |
| 16→32 | .813539 | 1.295799 | 1.113216 | -.257649 | .855567 |
| 32→64 | 1.101725 | 1.616906 | 1.335000 | -.222348 | 1.112652 |

Every pair is below 1.5 even when the endpoint/remainder cross term is omitted;
weakening cancellation lowers every observed order further. Subtracting E from
these fields would leave supplementary orders .7092,1.2958,1.6169: even that
hypothetical diagnostic change would not cure the first two failures.

The first row contains 69.44%, 5.58%, 18.01%, 50.62% of frozen squared norm.
In particular n=16 has exceptionally strong first-row cancellation. A first-row
norm fraction is therefore not a reliable measure of endpoint error magnitude.
A single O(h) row gives O(h^(3/2)) fixed-band RMS only asymptotically, assuming
smooth limiting derivatives and no leading cancellation. It does not imply
orders of 1.5 on these finite grids.

## What remains after endpoint isolation

To quantify the remainder, use a sequential supplementary split: matched radial
closure to seven-point radial first derivatives with the SAME second-order
angular derivatives, then to Fourier angular derivatives. Reconstruct physical
tensor divergence each time, on unchanged psi,W and original matter. This gives
F = E + R + A + Q, exactly up to roundoff. Q is the seven-point/Fourier physical
residual. R is the remaining radial change; A is the angular change. Mixed errors
depend on this order of replacement; these are not unique orthogonal components.
Full pairwise inner products are saved in budget.json.

| n | ||R|| | ||A|| | ||Q|| |
|---:|---:|---:|---:|
| 8 | 1.354315e-3 | 9.735076e-4 | 6.631803e-5 |
| 16 | 7.506655e-4 | 5.245477e-4 | 1.064979e-5 |
| 32 | 3.117249e-4 | 1.817358e-4 | 3.980951e-6 |
| 64 | 1.043906e-4 | 5.256472e-5 | 5.301479e-7 |

Both remaining radial and angular truncation are material, especially on coarse
grids. Q is much smaller than the frozen residual, consistent with the prior
supplementary accuracy result. This split and the exact endpoint identity
provide a quantitative mechanism for all three failures. They do not prove a
continuum impossibility theorem or give absolute solution error bounds.

## Improvement decision and preservation

Matching endpoint derivative errors is mathematically justified as a diagnostic
accuracy improvement, and the existing matched-closure and direct-Hessian
experiments already pursue consistent solver derivatives. Changing the frozen
audit is prohibited and would not address the measured coarse-grid remainder.
Adding a solver forcing designed to cancel E would deliberately introduce the
opposite O(h) first-row defect in an already much more accurate vector solution.
These measurements give no accuracy-based justification for that change.
No new solver modification or fresh solve is warranted by this analysis. This
is a concrete negative result, not a claim that no acceptable method can exist.

New analysis code and evidence are isolated on the research branch. All 105
original repository NPZ hashes (96 research plus 9 older reproduce archives)
match the prior inventory. No prior logs, failed results, checkpoints, histories,
source files, or frozen criteria were overwritten. Twenty polynomial checks
verify the cubic coefficient, quartic coefficient, and exact endpoint identity
on all four grids. Saved-field identities are checked by the analysis itself.

Reproduce from repository root:

    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/analyze_run43_endpoint_budget.py
    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/check_run43_endpoint_budget.py

These commands write only this study's budget/check JSON, not existing evidence.
Results: research/results/endpoint-error-budget-20261103/budget.json and
checks.json. Their SHA-256 values and the analysis source hashes are recorded in
SHA256SUMS in the same directory. No promotion, push, or merge was performed.
