# Run 43 independent angular refinement control

The angularly refined reference squared supplies **94.80%, 92.38% and 98.76%** of the three failed inner momentum gate deficits. The complete frozen hybrid audit reproduced **9/12** with bit-identical hybrid inputs.

## Method

For each frozen audit grid \(N=8,16,32,64\), the production solves used \(8N\) uniform radial nodal intervals, equivalent to \(2N\) continuous quartic elements, matching the finest archived reference for that audit grid. Each solve used \(2N\) or \(4N\) nodes in each angle. The independent radial Galerkin operator, physical sources, geometry, wall traces and angular stencil formulas come from the archived reference. In the equivariant reduction, refining \(\theta\) changes the centered derivative, face flux and damping, and refining \(\phi\) changes the rotation-generator derivative and the azimuthal damping factor. Restriction selects coincident nodes in all three directions. The scalar is constant in \(\phi\), and the cylindrical-frame displacement is rotated to Cartesian at the coincident audit azimuths.

The radial sensitivity check repeats each refined-angle solve with \(4N\) radial nodal intervals and compares the resulting residual vectors on the same audit nodes. The original-angle radial check uses the existing archive at both radial resolutions. The manufactured displacement check from the reference runs at each new solve grid. The matrix assembly and nonlinear iteration are reused verbatim. The sparse direct solves use a symmetrized sparsity graph for METIS nested dissection and SuperLU on the permuted matrix to limit fill.

The frozen order-two momentum auditor evaluates all sampled fields on the original audit grids. Write \(A=F_{4N}\), \(B=F_N-F_{4N}\) and \(C=H-F_N\), with \(F_N\) from the archived finest radial field. On each audit node, \(H=A+B+C\). Each reported squared norm and signed cross term uses the same hybrid weights \(r(R+r\cos\theta)\psi_{\mathrm{hybrid}}^6\), normalized over the indicated region. The first row is the first radial node inside the wall. The full inner band contains all interior audit nodes with \(r\leq0.10\).

For each term \(E\), the signed contribution to the \(p\geq1.5\) gate deficit is \(E_{\mathrm{fine}}-E_{\mathrm{coarse}}/8\). The six contributions add to the hybrid deficit. The first-row deficit diagnoses the row and the inner-band deficit is the frozen gate. The pair table includes changes in each contribution when using \(F_{2N}\), when substituting both half-density radial controls, and when substituting each half-density radial control independently.

## Gate attribution

The entries below are signed percentages of the hybrid gate deficit. \(AB\), \(AC\) and \(BC\) denote the signed doubled inner products. A negative share in a failed pair offsets part of the deficit. The first-row 8→16 deficit is negative, so its percentage signs are relative to a passing margin.

| Region | Pair | Hybrid order | Hybrid deficit | \(A^2\) | \(B^2\) | \(C^2\) | \(AB\) | \(AC\) | \(BC\) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| First row | 8→16 | 1.750 | −7.3002e−8 | 71.28% | 0.84% | −22.58% | 36.36% | 14.51% | −0.42% |
| First row | 16→32 | −0.616 | 3.9241e−7 | 79.98% | −0.039% | 0.16% | −2.00% | 21.75% | 0.15% |
| First row | 32→64 | −0.119 | 4.3670e−7 | 96.01% | −0.003% | −0.027% | −0.055% | 4.06% | 0.012% |
| Inner band | 8→16 | 0.475 | 5.0399e−7 | 94.80% | −0.13% | 0.11% | 0.66% | 3.99% | 0.57% |
| Inner band | 16→32 | 0.881 | 1.1276e−7 | 92.38% | −0.087% | −0.28% | −1.13% | 8.95% | 0.17% |
| Inner band | 32→64 | 1.123 | 1.6794e−8 | 98.76% | −0.046% | −0.17% | −0.95% | 2.37% | 0.034% |

On 8→16 the first row passes p≥1.5, while the full band fails and is controlled by \(A^2\). The leading correction in the band is the positive \(AC\) cross term. On 16→32, \(AC\) reinforces the failure by 8.95% in the band and 21.75% on the first row. The \(AB\) cross term offsets part of that contribution. On 32→64, \(A^2\) accounts for nearly the entire deficit in both regions. The positive \(AC\) and negative \(AB\) terms provide the largest corrections in the band. The six terms close against the hybrid squared norm to within 4.24e−22 across the eight decompositions.

## Refinement sensitivity

| \(N\) | Region | Weighted RMS of \(F_{4N}-F_{2N}\) | Signed RMS change | Radial change at \(2N\) | Radial change at \(4N\) |
| ---: | --- | ---: | ---: | ---: | ---: |
| 8 | First row | 3.575e−5 | −2.818e−5 | 1.747e−9 | 1.760e−9 |
| 8 | Inner band | 2.929e−5 | −2.300e−5 | 1.210e−9 | 1.219e−9 |
| 16 | First row | 1.124e−5 | −5.613e−6 | 3.738e−11 | 3.749e−11 |
| 16 | Inner band | 8.681e−6 | −4.622e−6 | 1.812e−11 | 1.817e−11 |
| 32 | First row | 3.058e−6 | +7.405e−7 | 4.276e−13 | 4.284e−13 |
| 32 | Inner band | 2.309e−6 | −7.847e−7 | 1.505e−13 | 1.515e−13 |
| 64 | First row | 8.221e−7 | +1.154e−7 | 8.731e−15 | 7.700e−15 |
| 64 | Inner band | 5.922e−7 | −1.316e−7 | 3.875e−14 | 3.021e−14 |

The signed change is \(\|F_{4N}\|-\|F_{2N}\|\). The radial columns use the norm of the vector change between the \(4N\) and \(8N\) radial meshes at fixed angles. They are at least four orders of magnitude smaller than the angular vector change. Halving both radial meshes shifts any band term by at most 0.000075 percentage points of the hybrid deficit. On the finest grid the radial differences provide an empirical sensitivity at floating-point scale.

Using \(F_{2N}\) in place of \(F_{4N}\) moves the band share of \(A^2\) by +0.118, −0.231 and −0.191 percentage points. That leaves the dominant attribution intact. The first-row 8→16 split between \(A^2\) and \(AB\) is more sensitive, moving by +7.43 and −7.08 percentage points of the passing margin. Some small band contributions have appreciable relative sensitivity. For example, the 32→64 \(AB\) share shifts from −0.946% to −0.772%, and \(B^2\) shifts from −0.0465% to −0.0297%. The measured angular vector change is about one fifth of \(\|F_N-F_{4N}\|\). These are finite-resolution comparisons. The remaining angular error at \(4N\) has no certified bound, and the residual of the \(4N\) field under the frozen auditor can contain both frozen differentiation error and remaining field error.

## Evidence and reproduction

- [Decomposition and refinement sensitivities](results/decomposition.csv) for each first row and inner band. The angular change is reported both as the RMS of the vector difference and as the signed change in RMS. The squared terms and crosses close against the independently computed hybrid energy.
- [Signed gate contributions and sensitivity substitutions](results/pairs.csv), [solver histories and generated field hashes](results/histories.json), [input and code hashes](results/provenance.json), and [all 12 frozen hybrid checks](results/hybrid-audit.json).

Run from the repository root with Python 3.12, NumPy 2.3.5, SciPy 1.17.0, PyMetis 2025.2.2 and mpmath 1.3.0. Use fresh output directories. The generated angular fields are scratch outputs. The analysis reads the existing radial archive in place. If the original hybrid archive is available, extract the hybrid directory and use it directly. The commands shown reproduce the four hybrid files byte for byte against the recorded hashes.

~~~bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pymetis==2025.2.2 mpmath==1.3.0
export OPENBLAS_NUM_THREADS=1
python research/run43_hybrid/physical.py /tmp/run43-hybrid
python research/run43_angular/solve.py /tmp/run43-angular
python research/run43_angular/analyze.py /tmp/run43-angular /tmp/run43-hybrid /tmp/run43-angular-results
~~~
