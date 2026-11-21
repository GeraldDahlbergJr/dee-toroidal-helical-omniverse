# Run 43 coupled momentum recovery experiment

**The new physical fields scored 9/12 in the original frozen suite.** The three inner-band momentum orders were 0.855455, 0.971827 and 1.127653, each below the original 1.5 threshold. All six Hamiltonian checks and all three outer-band momentum checks passed. The coupled solve converged at all four resolutions. The continuum evaluation showed finer-grid improvements in momentum residual, with a larger displacement error against the refined reference than the hybrid solver on every tested grid. [Unchanged frozen evaluation](results/fields/independent-audit.json) · [Numerical comparisons](results/comparison.csv)

## Discretization

Starting commit `84c859e1fda86d0dd35b737a0feddd3726aca7e1`. The experiment uses the hybrid study's azimuthally equivariant sector, toroidal geometry, matter projections, coupling constant and Hamiltonian face-flux stencil. The momentum equation is the divergence of the conformal strain with forcing $8\pi G\psi^6 S$. The Hamiltonian uses the squared strain from the same nodal first derivatives, including its existing matter energy and scalar angular stabilization. The physical Dirichlet data are $W=0$ and $\psi=1$ on both radial walls, with periodic angles.

For every ordered pair of outer and inner directions $(i,j)$, the separate sparse assembly adds $L(D_iD_jW,A_j)A_i + L(D_jW,\partial_i A_j)A_i$. Here $A_j$ are the Cartesian metric vectors and $L(v,a)=v\otimes a+a\otimes v-2(v\cdot a)I/3$. The interior pure radial term uses the three-point compact second derivative. The endpoint diagnostic rows use the reflected four-point second-order closure with weights $[2,-5,4,-1]/h^2$. The endpoint equations impose Dirichlet data by elimination. The adjacent interior rows use the compact stencil with the actual wall values. The periodic centered poloidal first derivative and the centered azimuthal spin factor are carried over from the diagnostic, as is composition for the pure angular second derivatives. The resulting operator is a strong collocation discretization. [Implementation](solver.py) · [Original diagnosis](../run43_frozen_momentum/REPORT.md)

## Manufactured validation

Two analytic vector fields excite radial, poloidal and azimuthal terms, with zero walls in one case and nonzero lifted walls in the other. Analytic forcing comes from the independent cylindrical-frame divergence. A separate discrete-forcing test recovers the manufactured nodal fields to within $2\times10^{-12}$ in max norm. The free rows agree with the original combined substitution applied to the manufactured fields up to roundoff. Both walls are reproduced exactly in the solves.

| Case | Inner wall | First row | Interior | Last row | Outer wall |
| --- | ---: | ---: | ---: | ---: | ---: |
| Homogeneous forcing defect, 64→128 | 1.907 | 1.827 | 2.004 | 2.073 | 2.296 |
| Lifted forcing defect, 64→128 | 1.992 | 1.985 | 1.973 | 1.951 | 1.958 |
| Homogeneous solved field error, 64→128 | exact | 2.986 | 2.017 | 3.000 | exact |
| Lifted solved field error, 64→128 | exact | 2.950 | 2.019 | 3.015 | exact |

The estimated inverse one-norm stays near $5.6\times10^{-3}$ on the tested meshes, while the condition estimate grows from about 50 to 3,579 over the physical grids. This is finite-grid evidence for bounded response to forcing in the tested sector. The composed poloidal derivative still has a Nyquist zero, and the test does not establish an angular ellipticity bound or stability outside the equivariant sector. [Manufactured tests](manufactured.py) · [Errors and solver diagnostics](results/manufactured.json)

## Coupled physical solves

Each grid converged in two vector–scalar updates after the initial momentum solve. The scalar residual is evaluated after refreshing the vector and strain. The last vector change was at most $2.84\times10^{-17}$.

| N | Momentum max residual | Momentum L2 residual | Hamiltonian max residual | Inner-band RMS $W$ difference from hybrid | Inner-band RMS $W$ error against 4N reference |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 4.16e−17 | 1.22e−16 | 6.35e−11 | 6.39e−7 | 7.53e−7 |
| 16 | 1.98e−16 | 1.10e−15 | 2.62e−12 | 1.36e−7 | 1.74e−7 |
| 32 | 1.62e−15 | 1.16e−14 | 1.15e−11 | 3.08e−8 | 4.14e−8 |
| 64 | 7.24e−15 | 1.18e−13 | 6.56e−11 | 7.38e−9 | 1.01e−8 |

Inner-band RMS differences in $\psi$ from the hybrid fields were $1.13\times10^{-12}$, $1.23\times10^{-12}$, $4.31\times10^{-13}$ and $1.18\times10^{-13}$. The corresponding errors against the refined reference were $3.89\times10^{-6}$, $9.92\times10^{-7}$, $2.46\times10^{-7}$ and $6.09\times10^{-8}$. Per-grid histories and field hashes accompany the raw arrays. [Fields and histories](results/fields) · [Regional field differences](results/comparison.csv)

## First row and full inner band

Orders below are $\log_2(E_N/E_{2N})$ for volume-weighted RMS momentum residuals. The first row means $r=R_{\min}+h$. The full inner band includes every interior row with $r\le 0.10$. The diagnostic 12/12 result applies the substituted momentum operator to the original hybrid fields. Its Hamiltonian checks come from the frozen evaluator.

| Fields and evaluator | First row 8→16 / 16→32 / 32→64 | Full inner band 8→16 / 16→32 / 32→64 | Score |
| --- | --- | --- | --- |
| Hybrid, frozen | 1.750 / −0.616 / −0.119 | 0.475 / 0.881 / 1.123 | 9/12 frozen |
| Hybrid, combined substitution | 1.770 / 2.033 / 2.090 | 1.868 / 1.568 / 1.630 | 12/12 diagnostic |
| Coupled, frozen | 1.955 / −0.313 / −0.144 | 0.855 / 0.972 / 1.128 | **9/12 frozen** |
| Hybrid, continuum | 1.856 / 2.753 / 0.597 | 1.507 / 1.620 / 1.749 | — |
| Coupled, continuum | 2.000 / 3.207 / 1.744 | 2.167 / 2.414 / 2.066 | — |
| Refined reference, frozen | 1.635 / −0.225 / −0.210 | 0.556 / 0.908 / 1.112 | — |

For the coarse hybrid and coupled arrays, each radial interval is reconstructed from five nearby coarse samples. The unchanged independent evaluator then differentiates each quartic trace at the coarse joins, differentiates the poloidal trigonometric interpolant and takes the exact azimuthal spin derivative. The refined reference uses its native quartic traces on 8N radial intervals and its 4N poloidal samples, restricted to the same evaluation nodes after differentiation. Each field uses its own $rq\psi^6$ weights and the original analytic source. The comparison table records sensitivity to shifting the radial interpolation stencil. At N=64 the shift changes the inner-band continuum residual vector by $4.44\times10^{-6}$ RMS for the coupled fields. The alternate reconstruction also gives lower coupled inner-band RMS at N=32 and 64, with a 33% reduction at N=64. The first row uses the same wall stencil under both biases, so that shift provides no independent sensitivity measurement there. [Evaluator](evaluate.py) · [Independent tensor differentiation](../run43_asymptotic/operator_eval.py)

The coupled continuum inner-band RMS is $7.31\times10^{-4}, 1.63\times10^{-4}, 3.06\times10^{-5}, 7.30\times10^{-6}$. It exceeds the hybrid residual at N=8 and 16 and falls below it at N=32 and 64. At N=64 the inner-band reduction is about 37% and the first-row reduction about 20%. The coupled displacement is still farther from the 4N reference than the hybrid displacement. These results distinguish differentiated residual convergence from field accuracy.

## Remaining obstruction

The same assembled momentum equation evaluated through the earlier combined diagnostic gives inner-band RMS at most $1.87\times10^{-15}$ on the new fields. The frozen implementation instead composes the radial first difference twice. At the first interior row its Taylor error includes $hW_{rrr}/4$, while its other interior rows have a wider second-order stencil. Its discrete differentiation of metric products also differs from the analytic product rule. Their cancellation varies with resolution. Removing just the radial discrepancy on the new N=8 fields *increases* the frozen-to-continuum inner-band gap from $8.63\times10^{-4}$ to $1.76\times10^{-3}$. At N=64 the radial discrepancy accounts for 94.8% of the inner-band gap's signed projection. Removing it reduces the first-row gap from $6.95\times10^{-4}$ to $2.74\times10^{-4}$ and the full-band gap from $2.00\times10^{-4}$ to $9.30\times10^{-5}$.

At N=64 the frozen first-row RMS is $7.03\times10^{-4}$, against $7.72\times10^{-6}$ for the independent continuum reconstruction. That first-row contribution enters a band with a growing number of rows as N increases. The measured band order therefore masks a first-row residual that is still increasing over the two finest grid pairs. Once a smooth first-order row defect dominates asymptotically, volume averaging over a fixed-width band supplies an additional half order as the first-row weight shrinks like $h$. The current grids give no evidence of meeting the finite-grid gate. The refined reference itself shows the same frozen first-row behavior and fails the three frozen inner-band momentum order thresholds. The obstruction is the resolution-dependent mismatch between frozen differentiation and the continuum-consistent discretization, with radial wall truncation dominant at fine resolution. Tightening the nonlinear solve would leave this operator discrepancy orders of magnitude above the solve residual.

## Reproduction and provenance

Run from the repository root with Python 3.12.14. The wrapper validates the manufactured problem before generating the physical fields, runs the untouched auditor in its own process, and replays the prior hybrid and 4N reference comparison with the earlier hash checks. Use a fresh output directory.

~~~bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pymetis==2025.2.2 mpmath==1.3.0
work=$(mktemp -d)
research/run43_coupled_recovery/reproduce.sh "$work"
diff -r research/run43_coupled_recovery/results/fields "$work/coupled"
diff research/run43_coupled_recovery/results/manufactured.json "$work/manufactured.json"
diff research/run43_coupled_recovery/results/comparison.csv "$work/comparison.csv"
~~~

The original auditor Git blob is `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`. The original field hashes and archive digests remain in the earlier provenance manifest. The replay reproduced the original diagnostic result directory byte for byte. [Original hashes and provenance](../run43_asymptotic/results/provenance.json) · [Reproduction wrapper](reproduce.sh)
