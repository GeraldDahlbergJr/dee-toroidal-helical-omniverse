# Run 43 targeted continuum auditor v2 validation

**Recommendation.** The v2 diagnostics are ready for prospective validation on resolved grids with the three rules below. A prospective acceptance gate still needs a declared physical residual tolerance and evidence that the independent evaluation is accurate at that tolerance. The original frozen result is 9/12 for both hybrid and coupled fields. This study evaluated manufactured inputs and archived data with no new physical solves.

## Roundoff and convergence

For each band, use the original volume-weighted RMS. Form the pointwise operator scale
$$
 B = \frac{8}{3}\psi^{-10}\|W\|_{\infty,2}
 \sum_{i,j}\left(b_{ij}|A_i||A_j|+b_j|A_i||\partial_i A_j|\right)
 + |8\pi G\psi^{-4}S|.
$$
Here $b_i$ and $b_{ij}$ are absolute stencil row sums or upper bounds on those sums. The radial first and second sums are $1/h_r,4/h_r^2$ inside and $4/h_r,12/h_r^2$ at the walls. The angular first sums are $1/h_\alpha$, and mixed and composed angular sums use products. The coefficient $8/3$ bounds the three terms of the vector contraction by the triangle inequality. This scaling tracks the amplification of perturbations in the **inputs** through the uncancelled operator, including when the output nearly vanishes.

Set the floor to $F=64\epsilon_{64}\|B\|_{\mathrm{RMS}}$. This is a conservative diagnostic margin. Float64 versus extended-precision evaluation on identical sampled geometry and source differed by at most $0.00172\epsilon_{64}\|B\|_{\mathrm{RMS}}$ in the manufactured calibration and $0.0156\epsilon_{64}\|B\|_{\mathrm{RMS}}$ in the archived array checks. An alternating one-ULP change to the sampled vector gave a residual change of at most $0.204\epsilon_{64}\|B\|_{\mathrm{RMS}}$ in the manufactured cases. The calibration covered 16–64 intervals and vector amplitudes from $2^{-20}$ to $2^{20}$. Exact discrete balance, an injected signal near $4F$, and a prescribed second-order signal landed in their expected classes. The low angular mode gave independently measured truncation-error orders of 1.979–1.995. [Measurements](measurements.json)

All 24 coupled interior region/grid measurements were roundoff-limited. The inner-band residual rose from $1.72\times10^{-17}$ to $1.87\times10^{-15}$, while the floor rose from $9.50\times10^{-15}$ to $5.52\times10^{-13}$. At 64, changing evaluation precision moved the inner-band residual by $3.65\times10^{-18}$ RMS. Most of the tiny residual therefore persists for the stored samples and the discrete equations. It is compatible with the rounding of those samples and the discrete solve. [Measurements](measurements.json)

## Angular sensitivity

The centered first derivative has response ratio $\sin(kh)/(kh)$. Its composition has squared response ratio. The independent angular FFT differentiates each pure second derivative directly with $-k^2$, including the even-grid Nyquist cosine. The real-grid first derivative has zero Nyquist coefficient under the cosine convention. The manufactured check compares the full operator response to an analytic Cartesian Hessian on each axis separately. [Test code](validate.py)

| Mode at 64 nodes | v1 response / analytic RMS, θ | v1 response / analytic RMS, φ |
| --- | ---: | ---: |
| Low, (k=1) | 0.9970 | 0.9971 |
| Intermediate, (k=16) | 0.4053 | 0.4054 |
| Near Nyquist, (k=31) | 0.000975 | 0.000913 |
| Nyquist, (k=32) | 0 | 0 |

The spectral relative error stayed below $10^{-11}$ on all 24 manufactured cosine tests. At exact Nyquist, a radially constant checkerboard has a nonzero continuum divergence but v1 returns zero at every node. With a quadratic radial envelope that vanishes at both walls, the v1 output equals the checkerboard sign times the radial-only response exactly. At 64 its weighted interior RMS is 56.9 in either angle, while the spectral-angular results with the same radial stencil are 12,752 in θ and 316 in φ. The fixed radial grid and growing angular frequency expose the loss of control of angular curvature by the v1 residual. The experiment measures residual sensitivity and rules out a uniform bound controlling angular curvature by v1 in this test family. Full coupled solver stability and iteration behavior are outside this result. A Nyquist sine vanishes on the original nodes and reaches unit amplitude on a half-cell shift, so FFT evaluation on the original samples cannot establish resolution of both phases. [Measurements](measurements.json)

## Prospective rules

1. **Roundoff.** Report $R$, $F$, and the norm and region used. Label $R\le F$ roundoff-limited and $F<R\le10F$ near-floor. Compute a v1 order only when both residuals exceed $10F$ on comparable grids. A floor-limited residual supports a statement about discrete balance. Recheck the floor on materially different operator implementations, precisions or scales.
2. **Angular resolution.** Screen native $W$ and $\psi$ samples along both angles. Compute the Parseval RMS fraction of $k^2\widehat f$ in modes with squared centered response below 0.90. Flag a fraction above 0.05 or a Nyquist amplitude fraction above $64\epsilon_{64}$. Ignore derivative signals at the FFT noise estimate $64\epsilon_{64}k_{\max}^2\|f\|_2$. These are proposed screening thresholds. A flag calls for an independent refined angular evaluation on a common comparison grid. A clear screen characterizes the sampled data. For a new field with possible aliasing, verify with a genuinely finer or shifted sample set that can expose the missing phase. Interpolating existing samples alone cannot supply this evidence.
3. **Continuum consistency.** Use an independent angular derivative and radial evaluation, with radial reconstruction or trace uncertainty checked by refinement. Use the same source, regions and reference weights for a given field comparison. Measure the **vector** gap between v1 and the independent result. For a convergence claim, require the independent residual to be above its own numerical floor and the refinement and reconstruction uncertainties to be small enough to resolve its decrease. For an absolute claim, predeclare a physical tolerance and verify that the independent residual plus its uncertainty fits inside it. Treat v1 alone as inconclusive if its gap exceeds the allowed audit error, even when v1 shows favorable orders or discrete balance.

## Archived application

| Field | N=64 inner-band v1 RMS | Independent RMS | Vector gap RMS | Reading |
| --- | ---: | ---: | ---: | --- |
| Hybrid | $1.29\times10^{-5}$ | $1.16\times10^{-5}$ | $7.08\times10^{-6}$ | Archived norms and independent decrease support a finite residual. The gap is 61% of the independent norm. |
| Coupled | $1.87\times10^{-15}$ | $7.30\times10^{-6}$ | $7.30\times10^{-6}$ | Discrete balance at the floor. The independent reconstruction gives last-pair orders 2.066 inner and 2.015 outer. |
| Reference N | $8.70\times10^{-6}$ | $3.16\times10^{-6}$ | $9.03\times10^{-6}$ | Above the floor on every tested region and grid. Gap is 2.86 times the independent norm here. |
| Reference 4N | $8.98\times10^{-6}$ | $1.99\times10^{-7}$ | $9.03\times10^{-6}$ | The gap is 45 times the archived independent norm. |

Coupled and native-N reference samples trigger the combined angular screen at 8, 16 and 32 and clear it at 64. At 32, $\psi$ in θ triggers both screens and $W$ also triggers the native-N reference screen, and the 64-node theta high-mode fractions are about $3\times10^{-5}$ for $W$ and $2.86\times10^{-4}$ for $\psi$. The archived 4N comparison shows about a 16-fold reduction of the N=64 inner independent residual relative to the native-N reference, while the v1 gap stays near $9.03\times10^{-6}$. Earlier radial controls document sensitivity of the independent residual to reconstruction. The solved fields cover the azimuthally equivariant sector. [Measurements](measurements.json) · [v1 trace comparison](../run43_continuum_auditor_v1/REPORT.md) · [Archived radial controls](../run43_asymptotic/REPORT.md)

Hybrid and 4N raw arrays were absent from the checkout and shared workspace. Their recorded v1 and independent RMS values are well above the floors measured on the available arrays. Their v2 input-based floors and spectral screens remain unverified, so the table uses archived comparisons and the machine-readable results leave those fields unset. The hybrid and coupled independent reconstructions use local quartics on coarse samples and have no independent radial error bound. The roundoff experiment holds sampled geometry and source fixed, and the screens and safety margin are empirical. [Archive hashes and measurements](measurements.json)

## Reproduce the targeted tests

Run from the repository root with NumPy 2.3.5 and a platform whose `np.longdouble` has greater precision than binary64. The recorded environment is in the JSON. The command verifies the starting-commit blobs, recorded file and member hashes and the frozen auditor before using archived results.

~~~bash
OPENBLAS_NUM_THREADS=1 python research/run43_continuum_auditor_v2/validate.py /tmp/run43-v2.json
diff -u research/run43_continuum_auditor_v2/measurements.json /tmp/run43-v2.json
~~~

If the original raw arrays are available, append `--hybrid /path/to/prior/hybrid --refined /path/to/prior/levels` to the first command. Each file must match its recorded digest. The starting commit is `59823a57cf1d83c5a792f874c79e891dc6793167`.
