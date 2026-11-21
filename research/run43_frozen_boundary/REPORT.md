# Run 43 frozen momentum boundary asymptotics

The unchanged coupled and hybrid audits each reproduced **9/12**. The coupled inner-band momentum orders were **0.855455, 0.971827 and 1.127653**, against 1.5. All six Hamiltonian checks and all three outer-band momentum checks passed in both reruns. The frozen radial composition has a first-order defect at the first interior row. For smooth fields with nonzero third radial derivative at the inner wall, it limits the auditor's continuum consistency in the weighted inner band to order 1.5. On the archived grids, the defect dominates the finest acceptance deficit, contributes to the middle deficit, and accounts for little of the coarsest deficit. Run 43 remains unrecovered under the original gate. [Coupled audit](results/coupled-acceptance.json) · [Hybrid audit](results/hybrid-acceptance.json)

## Radial composition and metric terms

Let $D_r$ be the order-two frozen radial first derivative on $r_i=a+ih$, with $a=0.05$, $b=0.25$, and $h=(b-a)/N$. The left wall uses $(-3,4,-1)/(2h)$ and the interior uses the centered difference. At the wall, composing takes $(-3D_0+4D_1-D_2)/(2h)$. At row one it takes $(D_2-D_0)/(2h)$. At rows $2$ through $N-2$ it applies two centered differences. The resulting second-derivative stencils and Taylor errors at the evaluation row are

| Row | Offsets | Coefficients times $h^2$ | $D_r^2 f-f''$ |
| --- | --- | --- | --- |
| Left wall | 0, 1, 2, 3 | $(5,-11,7,-1)/4$ | $3hf'''/4+5h^2 f^{(4)}/24+O(h^3)$ |
| First interior | −1, 0, 1, 2 | $(3,-5,1,1)/4$ | $hf'''/4+5h^2 f^{(4)}/24+O(h^3)$ |
| Rows $2$ to $N-2$ | −2, 0, 2 | $(1,-2,1)/4$ | $h^2 f^{(4)}/3+O(h^4)$ |

Reflection gives the right wall and right adjacent row, reversing the sign of the cubic term while preserving the quartic term. The orders follow by multiplying the coefficients by the powers of the offsets and dividing by the corresponding factorial. In the centered region, Taylor expansion through degree five gives the displayed $O(h^4)$ remainder for a field with bounded sixth derivatives. For comparison, the compact centered second difference used by the coupled solver has error $h^2f^{(4)}/12+O(h^4)$ on **all** free radial rows. [Frozen code](../../reproduce/diagnose_toroidal_boundary_convergence.py) · [Coupled code](../run43_coupled_recovery/solver.py)

Use Cartesian vector components and define $A_r=e_r$, $A_t=e_\theta/r$, $A_p=e_\phi/q$, $q=1+r\cos\theta$, and
$$
L(v,A)=v\otimes A+A\otimes v-\tfrac23(v\cdot A)I .
$$
For each ordered outer/inner pair $i,j$, the frozen term is $D_i L(D_j W,A_j)A_i$. Write $v=D_j W$. The principal term is $L(D_i v,A_j)A_i$. The metric term is the corresponding discrete product-rule remainder
$$
M_{ij}^h=\big[D_i L(v,A_j)-L(D_i v,A_j)\big]A_i .
$$
For a centered difference of spacing $s$, Taylor expansion and bilinearity cancel the terms involving $\partial_i v$ at leading order and $\partial_i^3 v$ at the next order, giving
$$
M_{ij}^h=\left[L(v,\partial_i A_j)+\frac{s^2}{6}
\{3L(\partial_i^2 v,\partial_i A_j)+3L(\partial_i v,\partial_i^2 A_j)+L(v,\partial_i^3 A_j)\}\right]A_i+O(s^4).
$$
For an outer radial wall row, replace $s^2/6$ by $-h^2/3$, with remainder $O(h^3)$. In comparing with the continuum include $L(v-\partial_j W,\partial_i A_j)A_i$, also second order for a smooth field. These expansions use smooth extensions in the outer direction. In $rt$, the inner angular derivative supplies a smooth radial function, including across the wall. In $tr$, the inner radial stencil is fixed on the evaluation row and gives a smooth angular function there even on the wall. In $rr$, $A_r$ is radially constant, so the discrete metric remainder vanishes exactly. The surviving radial error acts through $L(W_{rrr},e_r)e_r=(I+e_re_r^T/3)W_{rrr}$.

In particular $rt$ has $\partial_r A_t=-e_\theta/r^2$, and $tr$ has $\partial_\theta A_r=e_\theta$. Their continuum metric contributions are
$$
\frac{-e_\theta(e_r\cdot W_\theta)+(2/3)e_r(e_\theta\cdot W_\theta)}{r^2},\qquad
\frac{W_r+(1/3)e_\theta(e_\theta\cdot W_r)}{r}.
$$
They do not cancel identically. The two mixed **principal** contractions for $v=W_{r\theta}$ do cancel in part. Orthogonality gives
$$
L(v,A_t)A_r+L(v,A_r)A_t
=\tfrac13\{A_t(v\cdot A_r)+A_r(v\cdot A_t)\}.
$$
The same identity holds discretely with $v=D_rD_t W=D_tD_r W$. All ordered mixed pairs and the pure angular pairs have second-order truncation if both angular spacings are $O(h)$, with sufficiently smooth periodic fields and positive $r,q$. Thus no metric term can cancel the $O(h)$ coefficient in the asymptotic expansion. Finite-grid cancellation can be substantial because the higher-order coefficients involve inverse powers of $r$ and angular derivatives. The archived budgets resolve the $rt+tr$ and $tt$ metric terms by region. [Signed budgets](results/deficits.csv)

## Weighted-band limit

Work in exact arithmetic with a fixed smooth vector field and a smooth positive conformal factor, or grid samples of such fields. Assume sufficient bounded derivatives for the expansions above and angular quadrature convergence. The exact source is identical in the discrete and continuum evaluations, so it drops out of their difference. The weight is $w=rq\psi^6$, positive and bounded above and below on the shell. The band excludes the wall and contains $O(h^{-1})$ rows of a fixed physical band $a<r\le0.10$. Define the angular function
$$
B=\tfrac14\psi(a)^{-10}(I+e_re_r^T/3)W_{rrr}(a).
$$
The full momentum error on row one is $hB+O(h^2)$ and on the rest of the band it is $O(h^2)$ uniformly. Multiplying the audit sums by $h\Delta\theta\Delta\phi$ leaves their ratio unchanged. The row-one weighted numerator therefore contributes $h^3\int w(a)|B|^2\,d\theta\,d\phi+O(h^4)$, while the remaining rows contribute $O(h^4)$. The denominator tends to $J=\int_a^{0.10}\int w(r)\,d\theta\,d\phi\,dr>0$, with error $O(h)$ for these smooth weights and periodic trapezoid quadrature. Consequently
$$
E_h^2=C^2h^3+O(h^4),\qquad
C^2=J^{-1}\int w(a)|B|^2\,d\theta\,d\phi .
$$
If $C>0$, then $E_h=C h^{3/2}(1+O(h))$ and the measured doubling order tends to $3/2$. The row-normalized RMS is order one and the RMS of the remaining band is at least second order in the same setting. If $B=0$, the band is at least second order. The wall is outside the gate. Including it with a mass of order $h$ would add a second contribution of order $h^3$ to the numerator for a generic field.

This limit concerns auditor truncation on smooth fields. It is neither a universal cap on the residual of a mesh-dependent solver nor a guarantee of a finite-grid order of at least 1.5. To transfer the leading asymptotic to a family of numerical residuals, it is sufficient that the sampled continuum solution has nonzero $B$ and the additional frozen residual due to solver field errors is $o(h^{3/2})$ in the band. More generally one must control differentiated field errors and their possible cancellation. Nodal solution convergence alone does not provide this control.

## Manufactured fields

Sample $W=f(x)Q(\phi)(\cos 2\theta,\sin\theta,1+0.2\cos\theta)^T$ with $x=(r-a)/(b-a)$ and $Q$ the azimuthal rotation about the vertical axis. Both $f=\sin(\pi x)$ and $f=x(1-x)$ give homogeneous Dirichlet data and excite the angular and mixed terms. The former has nonzero wall third derivative and the latter has zero radial cubic defect. Use $\psi=1+0.05\cos\theta$ to exercise the conformal weighting and scaling. The manufactured error is the frozen operator minus independent analytic cylindrical-frame divergence, each scaled by $\psi^{-10}$. The reduced evaluation is equivalent by rotational symmetry and is checked against the full frozen 3-D operator at $N=8,16,32$. The table shows RMS errors, with separate normalization for each region. [All errors, leading terms and doubling orders](results/manufactured.csv)

| Field | N | Wall | First row | Remaining band | Full band |
| --- | ---: | ---: | ---: | ---: | ---: |
| sine | 8 | 376.62 | 69.077 | 90.032 | 81.73 |
| sine | 16 | 132.41 | 47.894 | 27.56 | 32.46 |
| sine | 32 | 47.218 | 23.169 | 7.3988 | 9.898 |
| sine | 64 | 18.892 | 9.1502 | 2.054 | 2.7666 |
| sine | 128 | 8.3554 | 3.5847 | 0.55682 | 0.7587 |
| sine | 256 | 3.9261 | 1.5088 | 0.14597 | 0.21201 |
| sine | 512 | 1.9037 | 0.68322 | 0.037427 | 0.061867 |
| sine | 1024 | 0.93753 | 0.32437 | 0.0094795 | 0.019058 |
| quadratic | 8 | 93.372 | 11.109 | 25.938 | 20.924 |
| quadratic | 16 | 29.09 | 8.4959 | 8.3236 | 8.3569 |
| quadratic | 32 | 7.7869 | 4.6363 | 2.2128 | 2.5268 |
| quadratic | 64 | 2.0095 | 1.6291 | 0.61431 | 0.68973 |
| quadratic | 128 | 0.51108 | 0.47749 | 0.16769 | 0.17987 |
| quadratic | 256 | 0.12896 | 0.12893 | 0.0442 | 0.045907 |
| quadratic | 512 | 0.032396 | 0.033476 | 0.011369 | 0.011594 |
| quadratic | 1024 | 0.0081192 | 0.0085279 | 0.0028845 | 0.0029134 |

The sine orders over 512→1024 are 1.022 at the wall, 1.075 on row one, 1.981 on the remaining band, and 1.699 on the full band. The last full-band order still has a bulk contribution. At $N=1024$ the leading-row-only contribution to band RMS is 0.015698 versus 0.019058 for the full error. Its 512→1024 order is 1.502, and the first-row weight fraction drops from 0.005214 to 0.002600. The first-row radial RMS is 0.307841 versus 0.307844 predicted by the cubic term. The quadratic full-band order is 1.993 over the last pair, with the pure radial difference at roundoff.

## Archived hybrid, coupled and reference fields

The coupled arrays are the committed physical outputs. The hybrid arrays were regenerated with matching archived hashes. The reference uses the archived $N$-angular fields and the digest-matched regenerated $4N$-angular fields, both with $8N$ radial intervals. Each case uses its own $rq\psi^6$ weights, original nodal data and analytic source. The reference is evaluated on native radial quartic traces and a native angular trigonometric interpolant before restriction. The hybrid and coupled continuum estimates use the same local quartic reconstruction as the coupled recovery report. They are estimates from sampled solutions and do not independently bound derivative reconstruction error. The frozen reduced residual is checked against the original full 3-D momentum evaluator for every case. The reference order comparisons apply the momentum threshold to the same bands and grids. [Regional measurements](results/archived.csv) · [Previous continuum sensitivity study](../run43_coupled_recovery/REPORT.md)

| Field | First row orders | Remaining band orders | Full band orders |
| --- | --- | --- | --- |
| Hybrid | 1.750 / -0.616 / -0.119 | -0.030 / 1.082 / 1.526 | 0.475 / 0.881 / 1.123 |
| Coupled | 1.955 / -0.313 / -0.144 | 0.309 / 1.157 / 1.553 | 0.855 / 0.972 / 1.128 |
| Reference N | 1.678 / -0.133 / -0.218 | 0.014 / 1.090 / 1.527 | 0.649 / 0.936 / 1.120 |
| Reference 4N | 1.635 / -0.225 / -0.210 | -0.106 / 1.070 / 1.520 | 0.556 / 0.908 / 1.112 |

### Signed attribution of each failed grid pair

Let (F) denote the frozen momentum residual on a given field and $E=\langle F,F\rangle_w$ its normalized weighted energy. An order at least 1.5 requires $D=E_f-E_c/8\le0$, where (f,c) denote fine and coarse. For an additive residual component (V), allocate $\langle V,F\rangle_w$ at each grid. This splits each cross term equally between its factors, and the grid allocations sum to (E) when the components sum to (F). The signed deficit contribution is $d_V=\langle V_f,F_f\rangle_w-\langle V_c,F_c\rangle_w/8$. Positive values contribute to failure and negative values offset it. Shares can exceed 100% under cancellation. The table gives $d_V$ and its share of (D) for the first-row pure radial defect (frozen minus reconstructed radial term, including its higher-order truncation) and for the entire first-row operator gap. These columns overlap. The last columns give the order after separately subtracting each measured defect from the residual on both grids, as a diagnostic on unchanged fields. Subtraction is nonlinear in energy and gives a different quantity from the signed allocation. [Complete budgets](results/deficits.csv)

| Field | Pair | Frozen deficit | First-row radial | Full row gap | Order after radial removal | Order after full row gap removal |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Hybrid | 8→16 | 5.040e-07 | +3.511e-08 (+7.0%) | -5.226e-08 (-10.4%) | 0.760 | -0.195 |
| Hybrid | 16→32 | 1.128e-07 | +7.858e-08 (+69.7%) | +3.363e-08 (+29.8%) | 1.283 | 0.998 |
| Hybrid | 32→64 | 1.679e-08 | +1.818e-08 (+108.3%) | +1.624e-08 (+96.7%) | 1.588 | 1.490 |
| Coupled | 8→16 | 4.417e-07 | +4.687e-08 (+10.6%) | -7.783e-08 (-17.6%) | 1.001 | 0.349 |
| Coupled | 16→32 | 1.009e-07 | +8.251e-08 (+81.8%) | +3.360e-08 (+33.3%) | 1.373 | 1.082 |
| Coupled | 32→64 | 1.641e-08 | +1.852e-08 (+112.8%) | +1.647e-08 (+100.4%) | 1.621 | 1.517 |
| Reference N | 8→16 | 4.805e-07 | -2.631e-08 (-5.5%) | -9.006e-08 (-18.7%) | 0.639 | -0.209 |
| Reference N | 16→32 | 1.028e-07 | +7.105e-08 (+69.1%) | +2.451e-08 (+23.8%) | 1.250 | 1.005 |
| Reference N | 32→64 | 1.642e-08 | +1.841e-08 (+112.1%) | +1.617e-08 (+98.5%) | 1.581 | 1.491 |
| Reference 4N | 8→16 | 4.778e-07 | -2.538e-08 (-5.3%) | -8.075e-08 (-16.9%) | 0.563 | -0.355 |
| Reference 4N | 16→32 | 1.042e-07 | +7.072e-08 (+67.9%) | +2.507e-08 (+24.1%) | 1.232 | 0.984 |
| Reference 4N | 32→64 | 1.659e-08 | +1.837e-08 (+110.8%) | +1.616e-08 (+97.4%) | 1.575 | 1.484 |

On the coupled 8→16 pair, the full row gap has a negative signed share, and its removal lowers the order. The mixed metric terms alone allocate −76.3% on row one and +95.5% on the remaining band, while the $tt$ metric allocations are +32.7% and −36.6%, respectively. This is finite-grid cancellation, consistent with the Taylor analysis. On 16→32 the first-row pure radial share is +81.8%, but the full row gap contributes only +33.3%. The remaining-band frozen order is itself below the threshold on both coarser pairs. On 32→64 the first-row pure radial allocation is +112.8%, and subtracting either radial or total first-row gap brings the coupled pair above 1.5. The hybrid and both reference variants also cross 1.5 on the finest pair after pure radial defect removal.

The $4N$ reference gives additional evidence independent of the coarse coupled solve. At $N=64$ its continuum band RMS is $1.99\times10^{-7}$, and its frozen band RMS is $2.00\times10^{-4}$. The frozen $4N$ reference still misses all three inner momentum orders. The coupled band at $N=64$ has continuum RMS $7.30\times10^{-6}$ and frozen RMS $2.02\times10^{-4}$. Both the reconstructed continuum residual and the hybrid–coupled differences show that physical discretization errors remain measurable, particularly on coarse grids. The frozen-versus-continuum gap on the refined reference and the manufactured proof identify the auditor as the primary obstruction. The evidence does not prove an upper bound on every possible solver's frozen score, establish a convergence theorem for the archived solver families, or attribute all three finite-grid failures to the first row alone.

## Next step

Develop a **separately versioned continuum-consistent momentum auditor** with independent manufactured qualification of the radial closure and metric differentiation. The evidence supports assessing that auditor on reference and physical fields before proposing another solver correction. Keep the original frozen score as the acceptance record. A new diagnostic cannot certify recovery under the original suite.

## Reproduction and integrity

Starting commit `a86b8f046908f852722d1d034f4e0711d5fc9466`. The reproduction uses Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, PyMetis 2025.2.2 and mpmath 1.3.0. Run from the repository root with a new scratch directory. The manufactured (1024) run uses the 2-D equivariant representation to avoid cubic storage.

~~~bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pymetis==2025.2.2 mpmath==1.3.0
work=$(mktemp -d)
research/run43_frozen_boundary/reproduce.sh "$work/replay"
diff -r research/run43_frozen_boundary/results "$work/replay/results"
~~~

The unmodified array audit checks the frozen auditor Git blob `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`, raw field hashes and boundary data. The coupled full audit is compared byte-for-byte with its archived copy. The hybrid full audit is checked against its prior SHA-256, and the entire earlier frozen momentum diagnostic output is recreated byte-for-byte. The script also checks the asymptotic manifest's source and output hashes, all named Run 43 reference archive array hashes and the regenerated $4N$ field hashes. It verifies the Run 42 ZIP digest, CRCs, and all 27 manifest member hashes. The branch diff from the starting commit contains only this research directory. [Verified digests](results/provenance.json) · [Replay script](reproduce.sh)
