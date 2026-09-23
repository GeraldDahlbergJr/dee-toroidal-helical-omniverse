# Time-dependent matter checkpoint

Run from the repository root:

```bash
python dee_radial_dynamics.py --config examples/dynamic_matter_config.json --output dynamic_matter_output
python -m unittest discover -s tests -v
```

The saved configuration evolves a regular radial Gaussian perturbation of all three fields in the manuscript's polynomial matter action on flat space. It uses cell-centered canonical momenta, finite-volume radial fluxes, reflecting outer flux, and fourth-order Runge–Kutta time stepping. The background is fixed and flat. The `rho`, `Theta`, and `Psi` amplitudes and their coupling are all active. The sample has **no angular winding** and does not prove finite-separation stability in general relativity.

## Numerical checks for the sample run

| Cells | Time step | Final time | Maximum relative energy drift | Minimum final `rho` |
| ---: | ---: | ---: | ---: | ---: |
| 120 | 0.005 | 1 | 4.30 × 10⁻¹² | 0.981317218 |
| 240 | 0.0025 | 1 | 1.70 × 10⁻¹³ | 0.981313331 |
| 480 | 0.00125 | 1 | 7.61 × 10⁻¹⁵ | 0.981312359 |

The RMS difference of the final `rho` profiles at coarse-grid cell centers is approximately `1.4611 × 10⁻⁶` for 120 versus 240 cells, and `3.6533 × 10⁻⁷` for 240 versus 480 cells, a ratio close to four. This is consistent with second-order spatial convergence for this sample. Interpolation and the finite domain are part of the comparison. The 240-cell run maintains a positive transport kinetic determinant, and the sampled radial and azimuthal null contractions are nonnegative. These checks address the coded matter equations and this chosen initial-value problem; they do not validate the physical DEE hypothesis.

The stress-energy profiles vary with time. Feeding one snapshot into the manuscript's *stationary* radial Einstein solver would silently change the problem, so this checkpoint stops at the matter sector. A gravitational comparison requires time-dependent tensor perturbation equations, a compatible gauge and constraint treatment, and documented boundary data. The archived weak-field CSVs report a different stationary construction with missing original source specifications; numerical agreement with them is neither claimed nor expected here.
