# Run 43 frozen inner-boundary consistency and asymptotic limit

The frozen momentum operator is consistent on the smooth toroidal shell. Its composed radial derivative is generically first order at the first interior row. On exact smooth samples, a nonzero third radial wall trace gives asymptotic order **1.5** in the original fixed-width weighted RMS when the angular grids are co-refined. The truncation law alone cannot guarantee the sign of the finite-grid margin over 1.5. It also does not bound the convergence order of every numerical solution sequence.

The coupled and hybrid fields retain their **9/12** frozen results. The finer coupled fields reduce the independently reconstructed continuum residual, while all three frozen inner-band momentum pairs fail. The first row controls the finest failure. It explains only part of the middle failure and actually offsets the coarsest failure. The physical grids have yet to exhibit the manufactured first-order row regime. [Derivation](DERIVATION.md) · [Physical measurements](results/physical.csv) · [Signed deficit allocation](results/deficits.csv)

## Manufactured evidence

The tests use $x=(r-0.05)/0.20$, radial factors $x-x^3$ and $x-x^4$, and the rotating Cartesian vector whose cylindrical components at azimuth zero are $(\cos(2\theta),\sin\theta,1+0.2\cos\theta)$. Both vectors vanish on the radial walls. The first radial factor gives a nonzero third derivative at the wall. The second gives zero third derivative there and an $O(h)$ third derivative at the first interior row. The positive factor $\psi=1+0.05x(1-x)\cos\theta$ gives the original $rq\psi^6$ weight a smooth nonconstant test. Analytic derivatives supply manufactured forcing independently of the difference stencils. Physical forcing is used for every physical comparison.

Measurements co-refine radial intervals and both angular counts from 8 through 2048. The equivariant evaluation uses the original radial and poloidal differences and exact rotation formulas for the frozen azimuthal differences. Direct 3D checks at 8 and 16 compare pointwise momentum at azimuth zero and the original full-band norm. At higher resolutions the calculation retains the band and two extra radial rows, which contain the full composed stencil of every reported row. No boundary closure at the end of that retained slice enters any measurement.

| Field | Quantity | 256→512 | 512→1024 | 1024→2048 |
| --- | --- | ---: | ---: | ---: |
| Cubic | Wall residual | 1.074 | 1.037 | 1.018 |
| Cubic | First-row residual | 1.223 | 1.119 | 1.060 |
| Cubic | First-row radial second derivative contribution, local RMS | 1.000 | 1.000 | 1.000 |
| Cubic | Radial second derivative contribution, full-band RMS | 1.504 | 1.502 | 1.501 |
| Cubic | Remaining rows, residual | 1.961 | 1.980 | 1.990 |
| Cubic | Full band, residual | 1.864 | 1.798 | 1.716 |
| Vanishing wall third derivative | First-row residual | 1.953 | 1.977 | 1.988 |
| Vanishing wall third derivative | First-row radial second derivative contribution, local RMS | 2.000 | 2.000 | 2.000 |
| Vanishing wall third derivative | Remaining rows, residual | 1.960 | 1.980 | 1.990 |
| Vanishing wall third derivative | Full band, residual | 1.986 | 1.993 | 1.997 |

The manufactured radial-only contribution includes the frozen tensor contraction and $\psi^{-10}$. For the cubic it vanishes analytically outside the first row within the band. The tiny measured radial-only values on those other rows are cancellation in floating point. At 2048, the first-row volume fraction is 0.001302, the cubic radial-only band RMS is 0.001069, and the full residual is 0.001328. The isolated contribution exhibits the predicted 1.5 order. The full cubic rate is still approaching it from above, with higher-order contributions evident even at 2048. The full quartic case approaches second order and demonstrates a positive finite-grid margin on an exact smooth manufactured field.

The first-row metric defect in the cubic test has order 1.993 on the finest pair. On that test at 2048 the first-row signed radial interaction with the two radial–poloidal metric terms is $+3.59\times10^{-5}$ in local mean-square units before doubling the cross term. On the coupled physical field at 64 the same interaction is $-1.98\times10^{-7}$, compared with radial-only mean square $8.16\times10^{-7}$. The signs and magnitudes depend on the field and resolution. [Manufactured measurements](results/manufactured.csv) · [Manufactured interactions](results/manufactured_interactions.csv) · [Physical interactions](results/physical_interactions.csv)

## Physical field comparison

The analysis uses the saved coupled fields, hash-matched hybrid fields, and the hash-matched 4N-angular, 8N-radial reference. Every residual compares values at the original audit nodes with the field's own conformal volume weights and the same analytic physical source. The continuum split uses the mean of one-sided quartic radial traces and Fourier angular differentiation. For hybrid and coupled samples, the quartics are local reconstructions from coarse nodal data. For the reference they use the native refined samples. The coarse reconstruction has no independent wall-jet accuracy bound. The prior coupled report's alternate stencil changes the full inner continuum residual vector by $4.44\times10^{-6}$ at 64 and shares the same first-row wall stencil. The reference comparison reduces dependence on coarse reconstruction, but is itself a numerical comparison. [Earlier reconstruction sensitivity](../run43_coupled_recovery/REPORT.md)

| Fields at N=64 | First row frozen RMS | First row continuum RMS | First row radial defect | First row metric defect | Inner frozen RMS | Inner continuum RMS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Hybrid | 6.989e−4 | 9.680e−6 | 8.991e−4 | 2.139e−4 | 2.031e−4 | 1.159e−5 |
| Coupled | 7.030e−4 | 7.723e−6 | 9.032e−4 | 2.127e−4 | 2.018e−4 | 7.301e−6 |
| Refined reference | 6.800e−4 | 2.903e−7 | 8.870e−4 | 2.116e−4 | 1.996e−4 | 1.993e−7 |

The coupled first-row frozen RMS grows on the two finer pairs, with orders −0.313 and −0.144. The first-row radial-only defect decreases with orders 0.201 and 0.510. The refined reference has similarly negative frozen row orders and a radial-only order of 0.442 on the finest pair. On that pair the remaining rows already give RMS orders of 1.526, 1.553 and 1.520 for hybrid, coupled and reference, respectively. The rising frozen first-row residual, declining radial defect and declining opposing metric contribution are compatible with a pre-asymptotic transition. They supply no measurement of the limiting physical wall trace or of the eventual physical convergence rate.

Define $E_N$ as the **squared** frozen full-band RMS. The deficit is $E_{2N}-E_N/8$. For each region, its contribution to $E_N$ is the local squared RMS times the fraction of the full-band weight in that region. This yields an exact, signed split across the first row and remaining rows. Term-level allocations divide each signed cross term equally between its two terms. The table's radial column is the first-row radial term's share under that convention. Negative shares offset failures. Shares above 100% indicate an offset by other rows or terms. The reference comparisons use the same momentum threshold as a diagnostic.

| Fields | Pair | Frozen order | Total deficit | First-row deficit | First-row share | First-row radial share |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Hybrid | 8→16 | 0.475 | 5.040e−7 | −7.298e−8 | −14.5% | 7.0% |
| Hybrid | 16→32 | 0.881 | 1.128e−7 | 3.304e−8 | 29.3% | 69.7% |
| Hybrid | 32→64 | 1.123 | 1.679e−8 | 1.651e−8 | 98.3% | 108.3% |
| Coupled | 8→16 | 0.855 | 4.417e−7 | −1.609e−7 | −36.4% | 10.6% |
| Coupled | 16→32 | 0.972 | 1.009e−7 | 3.012e−8 | 29.9% | 81.8% |
| Coupled | 32→64 | 1.128 | 1.641e−8 | 1.686e−8 | 102.8% | 112.8% |
| Refined reference | 8→16 | 0.556 | 4.778e−7 | −8.207e−8 | −17.2% | −5.3% |
| Refined reference | 16→32 | 0.908 | 1.042e−7 | 2.500e−8 | 24.0% | 67.9% |
| Refined reference | 32→64 | 1.112 | 1.659e−8 | 1.616e−8 | 97.4% | 110.8% |

## Conclusion and recommendation

The frozen operator is consistent for smooth fields. Its local first-row truncation error is generically first order, and the weighted band is exactly order 1.5 when the leading row coefficient is nonzero and the bulk is uniformly second order. A persistent positive asymptotic margin is impossible in that regime. A finite-grid margin is possible, as is a higher limiting order when the leading coefficient vanishes. There is no demonstrated positive margin for the physical fields on any of the three original frozen pairs.

The 9/12 results demonstrate a frozen discretization limitation for these fields. Even the refined reference exhibits the failure and a frozen RMS roughly one thousand times its continuum RMS at 64. The finest coupled failure is dominated by the first-row radial defect and its signed interactions. Finite-grid cancellation and the remaining rows matter on the coarser pairs. Differentiated numerical field error also remains, including sensitivity to continuum reconstruction. These data cannot exclude physical discretization errors or prove that an independently accurate solver can never pass the frozen checks. Physical observations are limited to the four available audit grids.

**Recommendation.** Introduce a separately versioned auditor with uniformly second-order consistency through the adjacent wall rows and a tested metric product rule, while preserving the frozen suite as the original benchmark. This targets the demonstrated evaluator error without requiring fields to compensate for the frozen stencil. The earlier substitutions motivate the design, but the new version would need its own verification. Run 43 remains at 9/12 under the original gate.

## Reproduction and provenance

Run from the repository root with Python 3.12.14. The wrapper regenerates hash-matched hybrid and refined-reference inputs through the existing replay, reruns both frozen acceptance audits, and verifies the coupled audit against its archive before measuring. Expect the physical gate to report three failed checks for each solver. Both previously archived arrays and the original evaluator stay in place.

```bash
python -m pip install numpy==2.3.5 scipy==1.17.0 pymetis==2025.2.2 mpmath==1.3.0
work=$(mktemp -d)
research/run43_inner_asymptotic/reproduce.sh "$work"
diff -r research/run43_inner_asymptotic/results "$work/results"
git show -s --format='%H %P' HEAD
```

The sole parent is `a86b8f046908f852722d1d034f4e0711d5fc9466` on `research/run43-frozen-momentum-20261121`. The isolated branch is `research/run43-inner-asymptotic-20261121`. Frozen evaluator Git blob `5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22`. Input digests are checked against the linked manifests. In the interaction tables, `signed_dot` is the local weighted dot product and `energy_band_normalized` applies the band mass fraction and doubles off-diagonal terms. The wall measurements are separate diagnostics and do not enter the gate. [Hybrid provenance](../run43_asymptotic/results/provenance.json) · [Refined reference output hashes](../run43_asymptotic/results/histories.json) · [Coupled output hashes](../run43_coupled_recovery/results/fields)
