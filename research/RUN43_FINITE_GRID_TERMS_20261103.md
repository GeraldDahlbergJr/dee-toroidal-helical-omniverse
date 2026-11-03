# Run 43: remaining finite-grid derivative terms

The completed endpoint budget commit `f110d62` was pushed to
`origin/research/run43-endpoint-error-budget-20261103`. This follow-up stays on
that dedicated research branch and uses only saved fourth-vector fields. No
baseline verification or solve was repeated.

The remaining radial error is not well represented by its leading smooth-field
term alone. A sizable quadratic derivative-composition term is concentrated at
the first interior row and cancels the radial linear term. Angular errors are
much better described by the leading centered-difference expansion by n=64.
These findings do not justify a solver change or candidate promotion.

## Exact factorial decomposition

Let G be the seven-point-radial/Fourier-angular reference gradient, R the
matched-radial-minus-reference gradient, and A the second-order-angular-minus-
Fourier gradient. These act on Cartesian components, including variable toroidal
basis coefficients. T(g,W) is the longitudinal tensor formed with gradient g;
D(g,T) is its Cartesian divergence. On fixed saved psi,W, the matched physical
residual minus the reference residual is exactly psi^-10 times:

    D(G,T(R,W)) + D(R,T(G,W))         radial linear
  + D(R,T(R,W))                       radial quadratic
  + D(G,T(A,W)) + D(A,T(G,W))         angular linear
  + D(A,T(A,W))                       angular quadratic
  + D(R,T(A,W)) + D(A,T(R,W))         radial/angular mixed.

Here “quadratic” refers to composing derivative errors, not nonlinear matter or
a source change. This explicitly exposes the mixed term hidden in the previous
sequential split. The original matter source cancels exactly in the differences.
All norms use the unchanged inner band and volume weights. Full Gram matrices
and row profiles are retained; components are not orthogonal.

| n | Radial linear RMS | Radial quadratic RMS | Angular linear RMS | Angular quadratic RMS | Mixed RMS |
|---:|---:|---:|---:|---:|---:|
| 8 | 1.856199e-3 | 7.927633e-4 | 9.402732e-4 | 1.717201e-4 | 1.411597e-4 |
| 16 | 1.002110e-3 | 4.019624e-4 | 5.196431e-4 | 1.223766e-5 | 1.781750e-5 |
| 32 | 3.819414e-4 | 1.324218e-4 | 1.814218e-4 | 7.815629e-7 | 1.527106e-6 |
| 64 | 1.176066e-4 | 3.227424e-5 | 5.254584e-5 | 4.890842e-8 | 1.096791e-7 |

The radial linear/quadratic cosines are -.8701,-.7764,-.6633,-.5265.
The quadratic norm is 42.7%,40.1%,34.7%,27.4% of the radial linear norm.
Its squared norm on the first interior row is 99.50%,98.55%,98.67%,98.72%.
Thus derivative composition near the closure remains material even after the
original endpoint term has been isolated. Calling this term uniformly O(h^4)
from smooth interior Taylor counting would be misleading on these grids:
derivative error operators also act on boundary-dependent error profiles.
This diagnosis is relative to the stated reference, not an absolute continuum
error measurement.

## Leading smooth-field predictions

Use R_lead(v)=er*h^2*partial_r^3(v)/6 with seven-point third derivatives, and
A_lead(v)=et*dt^2*partial_theta^3(v)/(6r)
          +ep*dp^2*partial_phi^3(v)/(6q), with Fourier angular derivatives.
Insert these perturbations in the LINEAR terms above, without fitted coefficients.

| n | Radial leading relative discrepancy | Angular leading relative discrepancy |
|---:|---:|---:|
| 8 | 20.79% | 16.36% |
| 16 | 19.83% | 3.757% |
| 32 | 17.43% | .905% |
| 64 | 13.92% | .223% |

Angular discrepancy decreases by roughly four per refinement, consistent with
the next centered-difference term. Radial finite-grid and closure effects remain
larger; the radial linear leading approximation additionally omits the sizable
quadratic term and mixed term listed above.

## Higher endpoint terms

At row 1, the exact forward third difference has the Taylor expansion

    delta3 W = W''' + (h/2) W'''' + (h^2/4) W'''''
               + (h^3/12) W'''''' + (h^4/40) W''''''' + ... .

Multiply by psi^-10*(h/4)*(I+er er^T/3) to obtain the endpoint contribution.
Using eight-point derivative estimates, including W'''' reduces its relative
error from 42.21%,30.29%,18.87%,10.60% to 21.18%,11.38%,4.56%,1.454%.
Adding W''''' reduces it to 4.50%,2.06%,.617%,.127%. Seven-point estimates give
comparable improvements (through W''''': 2.07%,1.05%,.367%,.090%). This supports
finite-h higher derivatives as the reason the leading endpoint term overpredicts
its amplitude. Completing the seven/eight-point polynomial series reproduces
the exact difference algebraically because all four samples belong to that
interpolant; that final closure is not independent continuum validation.

## Verification, limits, and preservation

The factorial identity closes with inner-band RMS discrepancies between
1.10e-17 and 1.31e-16 across the four grids. The initial numerical output and
subsequent profile-enriched output are both preserved separately:

- research/results/finite-grid-terms-20261103/terms.json
- research/results/finite-grid-terms-20261103-profile/terms.json

The initial file contains the same numerical decomposition before row profiles
were added. The final analysis script creates a new output directory and refuses
to overwrite existing results. Reproduce using a unique destination:

    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/analyze_run43_finite_grid_terms.py research/results/UNIQUE-finite-grid-output

These are supplementary derivative diagnostics on saved fields, with a finite-
order reference. They neither change nor replace the frozen 9/12 candidate
result. The accepted November 2 baseline remains untouched at its prior 12/12.
All 105 original archive hashes were compared with the existing inventory and
remain unchanged. Existing tracked solver/audit/source/checkpoint/evidence files
have no changes. New study files have SHA-256 hashes in the profile directory.

A follow-up method would need to improve whole-system boundary consistency and
accuracy while passing the same coarse-grid gate. This analysis identifies
substantial radial closure composition effects but does not establish a remedy
meeting both requirements. No source adjustment, fitted cancellation, new solve,
promotion, merge, or acceptance change was made.
