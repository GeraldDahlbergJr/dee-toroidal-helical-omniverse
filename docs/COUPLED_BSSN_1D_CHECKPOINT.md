# Symmetry-restricted BSSN evolution with nonlinear DEE matter

Run `PYTHONPATH=.:reproduce python reproduce/evolve_dee_bssn_variables_1d.py`
from the repository root. The numerical record is
`reproduce/coupled_dee_bssn_variables_1d_checkpoint.json`.

The evolved geometric variables are `phi = log(det(gamma))/12`,
`tilde gamma_ii = exp(-4phi) gamma_ii`, `K`,
`tilde A_ii = exp(-4phi)(K_ii-gamma_ii K/3)`, and the independent
contracted conformal connection `tilde Gamma^x`. The gauge is unit lapse,
zero shift, and periodic boundaries. The physical Ricci source is split
into conformal Ricci and conformal-factor Ricci; the former uses the
evolved connection in its derivative and Christoffel term. The BSSN trace
equation substitutes the Hamiltonian constraint, while the connection
equation substitutes the momentum constraint. The complete three-scalar
matter equations and independently checked source projections are coupled
at every RK4 stage.

| Points | Final Hamiltonian RMS | Final momentum RMS | Connection constraint RMS |
| ---: | ---: | ---: | ---: |
| 64 | 8.7648e-7 | 2.4212e-8 | 2.9226e-9 |
| 128 | 2.2050e-7 | 6.0714e-9 | 7.3289e-10 |
| 256 | 5.5211e-8 | 1.5190e-9 | 1.8336e-10 |

The adjacent-resolution maximum physical-state differences are 6.20245e-6
and 1.55247e-6, yielding order 1.99827 through `t=0.12`. The maximum
same-resolution difference from the separate ADM evolution is 6.79e-8,
1.72e-8, and 4.29e-9, respectively, also decreasing by approximately four
per refinement. `det(tilde gamma)-1` remains within 1.7e-15 and
`tr_tilde(tilde A)` within 9e-18. The evolved connection constraint is
independently compared with the derivative of the conformal metric.

The restrictions are substantial: all fields vary only in one Cartesian
direction; both geometric tensors are diagonal; the lapse and shift are
fixed; the box is periodic; and the run is short. The data are a small
sinusoidal transport-field perturbation near the finite radial branch.
They are not toroidal/helical initial data, and these results neither
establish long-time stability nor demonstrate a wormhole. Extending to
the intended configuration requires multidimensional geometry, an
appropriate gauge and boundaries, a constraint solve for that geometry,
and curvature and trapped-surface diagnostics.
