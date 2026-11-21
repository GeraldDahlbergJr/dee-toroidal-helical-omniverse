# Run 43 continuum momentum auditor v1

**The original frozen coupled result remains 9/12.** Its three inner-band momentum orders are 0.855455, 0.971827 and 1.127653 against 1.5. The new evaluator independently demonstrates second-order consistency on smooth fields at both walls, their adjacent rows, and the weighted bands. The hybrid arrays meet **12/12 new-auditor diagnostic checks**, combining the six original Hamiltonian checks and six v1 momentum checks. The coupled arrays satisfy the v1 stencil to roundoff, while their independently differentiated residual remains finite and converges. [Scores](results/scores.json) · [Frozen boundary study](../run43_frozen_boundary/REPORT.md) · [Coupled recovery](../run43_coupled_recovery/REPORT.md) · [Frozen momentum diagnosis](../run43_frozen_momentum/REPORT.md)

## Operator

Let $A_i$ be the Cartesian coordinate gradient in direction $i\in\{r,\theta,\phi\}$. The toroidal metric vectors are $e_r$, $e_\theta/r$, and $e_\phi/q$, with $q=R+r\cos\theta$. Write $L(v,a)=v\otimes a+a\otimes v-\frac23(v\cdot a)I$. The implemented divergence is
$$
\sum_{i,j} \left[L(D_{ij}W,A_j)A_i + L(D_jW,\partial_i A_j)A_i\right].
$$
Each mixed derivative is computed once and contracted in both orderings. The analytic metric table includes the azimuthal rotation of all three Cartesian basis vectors. For example, $\partial_\phi A_p=(\sin\theta\,e_\theta-\cos\theta\,e_r)/q$ and $\partial_r A_t=-e_\theta/r^2$. The implementation works on general 3-D Cartesian arrays. [Evaluator](auditor.py)

The pure radial term uses $(1,-2,1)/h^2$ on every interior row, with error $h^2W_{rrrr}/12+O(h^4)$. The wall uses $(2,-5,4,-1)/h^2$ and its reflection, each exact on cubics with leading error $-11h^2W_{rrrr}/12$. Mixed terms use tensor products of second-order first differences. The radial boundary first difference is quadratic one-sided and angular differences are periodic centered. For smooth inputs, the angular derivative acts on the smooth angular dependence of the radial truncation error, preserving second order at the walls. Pure angular derivatives compose the centered differences and are second order on smooth fields. Bounded smooth metric coefficients preserve these orders, and positive normalized weights preserve the second-order bound in the bands.

The residual retains $\psi^{-10}\operatorname{div}(LW)-8\pi G\psi^{-4}S$, using the frozen source and geometry. The norm uses the original $rq\psi^6$ weight and the inner and outer bands through $r=0.10$ and from $r=0.20$, respectively, excluding both endpoints. Endpoint compatibility is reported separately. The entry point reports momentum only. Hamiltonian values in the score are taken from the unmodified frozen audit.

## Independent manufactured qualification

The two tests refine all three directions at $N=8,16,32,64,128$, with $N$ radial intervals. Both exercise conformal scaling and volume weights with $\psi=1+0.03\cos\theta+0.02\sin\phi$ and use the original analytic source on both sides of the comparison.

- The homogeneous case uses a radial sine envelope and a nonconstant poloidal vector rotated about the vertical axis. Independent analytic cylindrical-frame derivatives yield the exact divergence. Both walls have zero displacement.
- The lifted case uses Cartesian plane waves with wavevectors $(1,2,3)$, $(-2,3,1)$ and $(3,1,-2)$ sampled at the toroidal Cartesian coordinates. This non-equivariant case varies with both angles and has nonzero displacement on both walls. The reference evaluates $\Delta W+\nabla(\nabla\cdot W)/3$ from the Cartesian Hessians, independently of the solver, finite differences and toroidal metric derivatives.

The radial stencil also differentiates polynomials of degree zero through three exactly at every row. For the last refinement, all manufactured RMS and maximum-norm orders exceed 1.8 across the eight measured regions. This manufactured qualification criterion leaves the historical acceptance threshold at 1.5. The table gives the 64→128 weighted RMS orders. [Qualification code](qualify.py) · [Manufactured error data](results/manufactured.csv)

| Case | Inner wall | Inner first | Inner remaining | Inner band | Outer wall | Outer first | Outer remaining | Outer band |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Homogeneous | 1.985 | 2.138 | 2.017 | 2.004 | 2.091 | 2.044 | 2.047 | 2.025 |
| Lifted Cartesian waves | 1.976 | 1.970 | 1.972 | 1.973 | 1.970 | 1.970 | 1.971 | 1.971 |

## Archived fields and independent continuum traces

Both discrete auditors use the raw 3-D hybrid and coupled arrays. Reference arrays have $8N$ radial intervals and either the archived $N$ or the digest-matched $4N$ poloidal samples. Both discrete auditors evaluate the reference fields after restriction to the audit grid. The independent continuum evaluator differentiates native radial quartics and Fourier interpolants before restriction, uses the exact azimuthal rotation generator, and averages the two radial traces at coarse joins. For the hybrid and coupled arrays, this evaluation uses local quartics reconstructed from coarse samples. Their estimates depend on reconstruction and carry no independent bound on derivative error. Each field uses its own original weights and source. [All regional RMS, maximum norms, orders and operator gaps](results/comparison.csv) · [Independent trace evaluator](../run43_asymptotic/operator_eval.py)

The table gives v1 inner-band RMS by region and pairwise observed orders. The coupled values at roundoff are included for reproducibility.

| Field | Inner region | N=8 | N=16 | N=32 | N=64 | Orders |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Hybrid | First row | 6.209e-04 | 1.821e-04 | 4.450e-05 | 1.045e-05 | 1.770 / 2.033 / 2.090 |
| Hybrid | Remaining band | 1.946e-04 | 9.726e-05 | 3.945e-05 | 1.300e-05 | 1.000 / 1.302 / 1.601 |
| Hybrid | Full band | 4.323e-04 | 1.184e-04 | 3.993e-05 | 1.290e-05 | 1.868 / 1.568 / 1.630 |
| Coupled | First row | 1.331e-17 | 4.369e-17 | 9.391e-17 | 2.341e-16 | -1.715 / -1.104 / -1.318 |
| Coupled | Remaining band | 1.955e-17 | 8.578e-17 | 3.779e-16 | 1.913e-15 | -2.133 / -2.139 / -2.340 |
| Coupled | Full band | 1.716e-17 | 7.944e-17 | 3.616e-16 | 1.872e-15 | -2.211 / -2.187 / -2.372 |
| Reference N | First row | 2.340e-04 | 1.065e-04 | 6.922e-05 | 2.679e-05 | 1.136 / 0.621 / 1.370 |
| Reference N | Remaining band | 3.200e-04 | 8.554e-05 | 2.203e-05 | 6.821e-06 | 1.903 / 1.957 / 1.691 |
| Reference N | Full band | 2.863e-04 | 8.995e-05 | 2.954e-05 | 8.695e-06 | 1.670 / 1.606 / 1.764 |
| Reference 4N | First row | 3.050e-04 | 6.844e-05 | 6.282e-05 | 2.578e-05 | 2.156 / 0.124 / 1.285 |
| Reference 4N | Remaining band | 4.144e-04 | 1.097e-04 | 2.617e-05 | 7.362e-06 | 1.918 / 2.067 / 1.830 |
| Reference 4N | Full band | 3.715e-04 | 1.030e-04 | 3.128e-05 | 8.981e-06 | 1.850 / 1.720 / 1.800 |

On the $4N$ reference the inner-band vector RMS gap between v1 and the independent continuum evaluator falls from $3.83\times10^{-4}$ to $9.03\times10^{-6}$ with orders 1.876, 1.730 and 1.802. The frozen gap on the same fields falls from $1.18\times10^{-3}$ to $2.00\times10^{-4}$. At $N=64$ the v1 gap is about 22 times smaller than the frozen gap. The independent reference residual itself is $1.99\times10^{-7}$ at that resolution. The v1 consistency error remains measurable at the tested resolutions and matters when interpreting the absolute discrete residual.

## Frozen acceptance and v1 diagnostic

| Fields | Original frozen acceptance | Original Hamiltonian checks | v1 momentum diagnostic | Combined new-auditor diagnostic |
| --- | ---: | ---: | ---: | ---: |
| Hybrid | 9/12 | 6/6 | 6/6 | **12/12** |
| Coupled | 9/12 | 6/6 | 0/6 | **6/12**, roundoff-scale threshold arithmetic |

For the coupled data, v1 uses the same interior momentum discretization as the solve. Its residual ranges from $1.72\times10^{-17}$ to $1.87\times10^{-15}$ in the inner band. Ratios at that scale are roundoff-sensitive. The independently reconstructed coupled inner-band residual instead decreases from $7.31\times10^{-4}$ to $7.30\times10^{-6}$ with observed orders 2.167, 2.414 and 2.066. Reference residuals are compared by operator and region without assigning a combined Hamiltonian score.

## Qualification scope and recommendation

**Use v1 as a versioned momentum diagnostic for future Run 43 comparisons on smooth, resolved toroidal grids.** Its wall and band consistency has independent analytic and refined-reference support. For solvers that share the stencil, include an independent continuum trace evaluation or another independently refined residual. The composed angular differences retain a Nyquist null mode, so smooth-field consistency establishes neither high-frequency sensitivity nor stability. The archived solved fields are azimuthally equivariant, while the non-equivariant validation is manufactured. An official acceptance gate would require a prospective protocol addressing underresolved modes and roundoff-scale residuals.

## Reproduction and preservation

The isolated branch starts at 3715241456898277e290fe135222582d538c182e. The replay checks that tracked files outside the v1 directory match that commit, reruns the original hybrid and coupled audits in separate processes, and compares the inherited frozen-boundary evidence byte for byte. It verifies the frozen auditor Git blob 5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22, coupled and regenerated hybrid array hashes, reference archive and member hashes, regenerated $4N$ arrays, earlier audits and provenance, and the Run 42 baseline archive and manifest member hashes. [Verified hashes and versions](results/provenance.json)

Run from the repository root with Python 3.12.14.

~~~bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pymetis==2025.2.2 mpmath==1.3.0
work=$(mktemp -d)
research/run43_continuum_auditor_v1/reproduce.sh "$work/replay"
diff -r research/run43_continuum_auditor_v1/results "$work/replay/results"
~~~

For a single input array file:

~~~bash
python research/run43_continuum_auditor_v1/auditor.py \
  research/run43_coupled_recovery/results/fields/fields_64.npz
~~~
