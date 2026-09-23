# Static radial sector: a constraint on the next DEE calculation

This note analyzes a **restricted** version of the manuscript's matter action, not the full toroidal–helical theory. It assumes a flat cylindrical background, smooth real fields `rho(r)`, `Theta(r)`, `Psi(r)`, no dependence on time, angle, or axial position, a regular axis, finite energy per unit axial length, vacuum behavior at infinity, and a positive-definite transport kinetic matrix. It takes the polynomial nonnegative potential in manuscript Eq. (14), with `lambda_V>0` and `v>0`.

## Transport fields

The transport equations (manuscript Eqs. (5)–(6)) become

```text
(r J_Theta)' = (r J_Psi)' = 0,
J_Theta = K_Theta Theta' + Lambda Psi',
J_Psi   = Lambda Theta' + K_Psi Psi'.
```

Regularity at `r=0` sets both integration constants to zero. An invertible positive-definite kinetic matrix then implies `Theta'=Psi'=0`. In this restricted sector, there is no smooth radial transport current.

## Radial field

The remaining static energy per unit axial length is

```text
E = 2 pi integral_0^infinity r [ (rho')^2/2 + V(rho) ] dr,
V(rho) = lambda_V (rho^2-v^2)^2/4 >= 0.
```

For a finite-energy stationary field, consider the allowed scaling `rho_s(r)=rho(s r)` with vacuum boundary conditions at infinity. Its gradient energy is independent of `s` in two spatial dimensions, while its potential energy scales as `s^-2`. Stationarity at `s=1` requires the integral of `V` to vanish. Since `V>=0`, a smooth continuous field has `rho^2=v^2` throughout; the vacuum choice at infinity gives `rho=v` everywhere. The source tensor of this vacuum is zero. This is the usual scaling obstruction for a nontrivial localized static scalar lump in this restricted two-dimensional sector.

Consequently, the positive, localized Gaussian example in `matter_profile_check.py` cannot be adjusted into a nontrivial regular on-shell profile merely by tuning its amplitudes. The numerical residuals demonstrate its failure; the argument above explains why a whole class of static radial attempts will fail under the stated assumptions.

## What changes the problem

A nonzero, localized DEE source requires changing one or more assumptions and solving the resulting equations with explicit boundary conditions. Possibilities to investigate include time or axial dependence, angular winding with a globally regular core and suitable field manifold/couplings, an imposed finite-domain boundary or external source, or a different potential or matter sector. These are distinct models and cannot be substituted into the current stationary radial Einstein calculation without checking their stress-energy conservation and the metric/gauge constraints. In particular, the manuscript's local `m phi` phase is not automatically regular at a global axis when its kinetic coefficient remains nonzero there.

**Review implication:** the existing weak-field CSVs remain archived exploratory results. The vacuum is an exact but trivial source check; it does not reproduce those CSVs. A new physical checkpoint needs an on-shell matter solution and a gravitational approximation compatible with its symmetries.
