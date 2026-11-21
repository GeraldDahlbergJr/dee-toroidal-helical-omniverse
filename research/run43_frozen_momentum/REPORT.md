# Run 43 frozen momentum operator diagnosis

The unchanged hybrid acceptance suite reproduced **9/12**, with the three inner-band momentum orders at **0.474714, 0.881283, and 1.123133** against 1.5. The dominant fine-grid discrepancy is the radial second derivative in the Cartesian vector Laplacian. At the coarsest grid the radial–poloidal metric terms are substantial and cancel other errors. A diagnostic that combines a compact radial second derivative with analytically differentiated metric coefficients meets **12/12 diagnostic thresholds** on the same hybrid fields. This score applies to the substituted diagnostic operator. The frozen suite result is 9/12.

Starting commit `5996f9f5de797aba7d67fc0025b7007ec40399fb`. The frozen auditor Git blob is `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`. The full independently regenerated frozen audit has SHA-256 `6c10d04f1ffbb9006b9544f8c2bb5669218a614f99d5e608db0d4fd6c4207de1`.

## Radial differentiation

The order-two frozen path applies a three-point radial first derivative to each Cartesian component of the displacement and then to each component of its conformal strain. The wall first derivative is the quadratic one-sided formula and the interior first derivative is centered. Each angular derivative is a periodic centered first difference. The gradient uses the Cartesian vectors $A_r=e_r$, $A_t=e_\theta/r$ and $A_p=e_\phi/q$, with $q=R+r\cos\theta$. The outer gradient contracts the second index of the symmetric trace-free Cartesian strain.

The frozen path implements its pure radial second derivative by composition. It has the following dimensionless rows after multiplying by $h^2$. The right wall reflects these rows.

| Row from inner wall | Offsets | Composed coefficients | Leading error |
| ---: | --- | --- | --- |
| 0 | 0, 1, 2, 3 | $[5,-11,7,-1]/4$ | $3hf'''/4$ |
| 1 | −1, 0, 1, 2 | $[3,-5,1,1]/4$ | $hf'''/4$ |
| 2–4 | −2, 0, 2 | $[1,-2,1]/4$ | $h^2f^{(4)}/3$ |

The compact centered second derivative on rows 1–4 uses offsets −1, 0, 1 and coefficients $[1,-2,1]$, with leading error $h^2f^{(4)}/12$. At the wall its one-sided three-point formula is first order. The six-point construction is fourth order at the wall and first interior row. Every row reproduces derivatives through degree two. The dimensionless manufactured tests use $x^3$, $e^x$ and $\sin(2x)$ on $[0,1]$ through 256 intervals. The cubic exposes the exact first-order defects $4.5h$ and $1.5h$ in the composed wall and first interior rows, respectively. Cubic errors in exact rows are floating-point noise. [Stencil coefficients](results/stencils.csv) · [Manufactured errors and orders](results/manufactured.csv)

| $e^x$ error at 64 intervals | Composed | Compact | Six-point |
| --- | ---: | ---: | ---: |
| Wall | 1.177e−2 | 1.577e−2 | 4.686e−8 |
| First interior | 4.020e−3 | 2.067e−5 | 4.453e−9 |
| Second interior | 8.397e−5 | 2.099e−5 | 6.840e−10 |
| First interior order, 32→64 | 1.041 | 2.023 | 4.049 |
| Second interior order, 32→64 | 2.045 | 2.045 | 4.044 |

## Term substitutions

For each ordered outer/inner direction pair $(i,j)$, the principal term is $L(D_iD_jW,A_j)A_i$, where $L(v,a)=v\otimes a+a\otimes v-2(v\cdot a)I/3$. The frozen metric term is the discrete remainder $[D_iL(D_jW,A_j)-L(D_iD_jW,A_j)]A_i$. In the continuum it is $L(\partial_jW,\partial_i A_j)A_i$. This split preserves both mixed derivative orderings and the tensor contraction. The exact radial–radial metric term vanishes because $e_r$ is independent of radius. The radial principal term is $W_{rr}+e_r(e_r\cdot W_{rr})/3$.

The independent evaluator uses the same toroidal geometry in a cylindrical orthonormal frame at $\phi=0$. It averages analytic quartic left and right strong traces at radial element joins, differentiates the poloidal trigonometric interpolant, and uses the exact equivariant rotation generator $J$ on vectors and $JL-LJ$ on tensors. The diagnostic sums match both original implementations and reproduce the prior asymptotic residual and gap values at the tested levels. Each paired evaluation uses identical restricted reference nodes for the vector, conformal factor, and analytic source. The normalized weight is $rq\psi_{\rm hybrid}^6$ at the same audit nodes. Equivariance makes the integrands independent of azimuth, so the sums can be evaluated at $\phi=0$. The frozen identity is also checked against the full three-dimensional auditor.

The sample includes the verified 1N-angular archive and the same 4N-angular fields used in the prior study, regenerated and checked against its field hashes. Both use 8N radial intervals. The 4N inner-band gap agrees with the prior finest measured gap to about 0.022% or better across the four grids. The table shows 4N results. A signed projection is $\langle\Delta T,\Delta M\rangle/\|\Delta M\|^2$. Joint and individual terms can exceed 100% under cancellation. “After exact rr” replaces only that term by its independently differentiated value.

| N | Frozen−continuum RMS | rr RMS | rr projection | rt + tr metric projection | tt metric projection | After exact rr |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1.176e−3 | 1.051e−3 | 14.5% | 248.2% | −114.1% | 1.444e−3 |
| 16 | 8.071e−4 | 7.871e−4 | 57.9% | 97.3% | −40.1% | 7.185e−4 |
| 32 | 4.310e−4 | 4.431e−4 | 86.8% | 32.2% | −13.8% | 2.440e−4 |
| 64 | 1.995e−4 | 2.044e−4 | 96.2% | 10.9% | −5.1% | 7.048e−5 |

At N=64 the radial vector Laplacian contributes 96.24% of the gap projection and the radial gradient-of-divergence term contributes 0.006%. Substituting the exact rr term reduces the first-row gap from $6.801\times10^{-4}$ to $2.132\times10^{-4}$. At N=8 the two radial–poloidal metric terms supply nearly all of the mixed-term projection. The combined mixed principal terms are small. The entire outer azimuthal derivative contributes at most 0.76% in signed magnitude across the inner-band grids. The term table records individual and grouped RMS values, observed orders, and the remaining gap and residual after each exact term substitution. Grouped rows overlap.

A separate signed allocation assigns each term $V$ its share $\langle V,H\rangle$ of the hybrid squared residual, with $H$ the frozen hybrid residual. This distributes each cross term equally between its two factors. Applying the original $E_{\rm fine}-E_{\rm coarse}/8$ deficit to those allocations gives the following shares. The other terms and cross-term effects are included in the CSV. The first interior row supplies most of the rr contribution to the two finer deficits. The other radial rows supply nearly all of the rr contribution to the coarsest deficit. [Term substitutions and orders](results/terms.csv) · [Signed deficit allocations](results/deficits.csv)

| Deficit pair | rr first row | rr other rows | rt + tr metric | tt metric |
| --- | ---: | ---: | ---: | ---: |
| 8→16 | 0.01% | 75.95% | 38.43% | −13.77% |
| 16→32 | 64.94% | 38.41% | −15.04% | 5.81% |
| 32→64 | 109.02% | −1.60% | −17.45% | 6.99% |

## Candidate diagnostics on unchanged fields

The compact and fourth-order rr variants replace only the pure radial principal term. The wall variants replace the first and last two rows of that term, including the endpoints. The metric product-rule variant uses the analytic derivatives of $A_j$ with the original frozen nodal first derivatives of $W$. The combined variant adds the compact rr correction and the metric product-rule correction. All retain the original physical source, nodal fields, boundary data, weights, Hamiltonian evaluation, and threshold.

| N | Frozen reference gap | Wall-only compact | All-row compact rr | Combined compact rr + metric |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 1.176e−3 | 2.158e−3 | 2.085e−3 | 3.835e−4 |
| 16 | 8.071e−4 | 1.215e−3 | 9.973e−4 | 1.044e−4 |
| 32 | 4.310e−4 | 4.755e−4 | 3.364e−4 | 3.149e−5 |
| 64 | 1.995e−4 | 1.532e−4 | 9.701e−5 | 9.029e−6 |

A closure-only change improves the finest gap and worsens the first three. The combined correction reduces the gap by 67.4%, 87.1%, 92.7% and 95.5%. On the hybrid fields it lowers the inner-band momentum RMS from $1.133\times10^{-3}$ to $4.323\times10^{-4}$ at N=8 and from $2.031\times10^{-4}$ to $1.290\times10^{-5}$ at N=64. [Reference residuals and gaps](results/reference.csv) · [Hybrid residuals and orders](results/hybrid.csv)

| Operator | Inner 8→16 | Inner 16→32 | Inner 32→64 | Checks passed |
| --- | ---: | ---: | ---: | ---: |
| Frozen acceptance suite | 0.474714 | 0.881283 | 1.123133 | 9/12 |
| Compact rr, wall rows only | 0.796892 | 1.315158 | 1.612326 | 10/12 diagnostic |
| Fourth-order rr, wall rows only | 0.720185 | 1.259490 | 1.578669 | 10/12 diagnostic |
| Compact rr, all rows | 1.022420 | 1.522873 | 1.761333 | 11/12 diagnostic |
| Fourth-order rr, all rows | 1.089887 | 1.583203 | 1.795058 | 11/12 diagnostic |
| Metric product rule | 0.719004 | 0.932674 | 1.164472 | 9/12 diagnostic |
| Compact rr + metric product rule | 1.868352 | 1.568025 | 1.629676 | 12/12 diagnostic |

The combined variant's smallest margin on the inner band is 0.068. The next experiment should run a coupled manufactured and physical solve with a separately assembled, wall-consistent momentum discretization informed by the combined diagnostic. Measure differentiated solution error against the quartic/Fourier reference and run the unchanged frozen acceptance suite on the new fields. The reference fields themselves have frozen inner-band orders of 0.556, 0.908 and 1.112 at 4N, so simply improving continuum field accuracy does not establish a route to frozen acceptance. A controlled companion solve that incorporates the audited operator would test whether the frozen orders can be recovered while also reducing the independently measured continuum residual. [Frozen checks](results/frozen-acceptance.json) · [Diagnostic checks](results/acceptance.json)

## Reproduction

Run from the repository root with Python 3.12.14. The wrapper reuses existing hybrid and 4N fields in the supplied work directory, runs missing solves serially, and independently reruns the unchanged acceptance runner. The diagnostic verifies the archive, auditor, source and field hashes against prior manifests and checks its reproduced residuals against the prior result table.

~~~bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pymetis==2025.2.2 mpmath==1.3.0
work=$(mktemp -d)
research/run43_frozen_momentum/reproduce.sh "$work" "$work/replay"
diff -r research/run43_frozen_momentum/results "$work/replay"
~~~
