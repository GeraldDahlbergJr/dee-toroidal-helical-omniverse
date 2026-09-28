# Three-dimensional helical ring: initial constraints and first evolution trial

Run `PYTHONPATH=. python reproduce/toroidal_helical_initial_data.py` and
`PYTHONPATH=.:reproduce python reproduce/evolve_dee_ring_3d.py`.
For the second hierarchy, run
`PYTHONPATH=.:reproduce python reproduce/evolve_dee_ring_3d.py --resolutions 16 32 64 --output coupled_dee_ring_3d_secondary_checkpoint.json`.

This trial uses the archived three-scalar nonlinear action. It **chooses**
a smooth finite-amplitude ring seed in a Cartesian periodic box of side 8:
`rho=1-0.025 F`, `Theta=0.035 F Re[((x+iy)/R)^m exp(2 pi i z/L)]`, and
`Psi=0.028 F Im[exp(0.4i)((x+iy)/R)^m exp(2 pi i z/L)]`, with `R=1.7`,
`m=2`, and a smooth localized envelope `F`. The integer angular mode is a
Cartesian polynomial and has no branch cut on the axis. These are regular
oscillatory real scalars, **not** globally wound phase fields; no toroidal
coordinate singularity or wormhole topology is imposed.

With zero scalar momenta, `S_i=0`, constant-mean-curvature isotropic
`K_ij=gamma_ij K/3`, and `gamma_ij=psi^4 delta_ij`, a periodic conformal
Hamiltonian solve adjusts `psi` and `K`. The discrete conformal equation
residual is below 6e-13 for 16/32/64 grids, with the momentum source zero.
An independent Ricci calculation on the resulting full tensor metric gives
initial Hamiltonian RMS residuals 2.89e-5, 1.63e-5, and 8.65e-6 on the
12/24/48 evolution grids. These are truncation residuals of a *different*
geometric discretization; the near-zero conformal residual must not be
presented as near-zero independent geometric residual.

The trial advances all six spatial metric components, all six curvature
components, and three coupled scalars with the 3-D ADM equations, unit
lapse, zero shift, periodic boundaries, centered differences, and RK4 to
`t=0.12`. The off-diagonal `K_xy` and `gamma_xy` become nonzero.

| Points | Steps | Final H RMS | Final momentum RMS |
| ---: | ---: | ---: | ---: |
| 12 | 1 | 2.91e-5 | 3.11e-6 |
| 24 | 2 | 1.42e-5 | 3.99e-6 |
| 48 | 4 | 4.61e-6 | 2.64e-6 |

The physical-state differences between adjacent grids are 0.01645 and
0.01278, an apparent order of just 0.364. The Hamiltonian residual falls,
but the momentum residual rises from 12 to 24. **This three-resolution
evolution fails the convergence gate.** The coarse grids poorly resolve
the ring and the number of time steps is small; those observations are
possible contributors, not a verified explanation. The 12/24/48 result
is retained as recorded. The geometry is a periodic ring of helical
matter, not a derived toroidal spacetime or a traversable wormhole.

The next calculation should independently check the full 3-D Ricci and
momentum operators against manufactured data, then repeat the same ring
at a demonstrably resolved hierarchy and with enough time steps. A
successful target-geometry run would additionally need an appropriate
gauge and outer boundary, a toroidal/helical spacetime seed and its
constraint solve, curvature invariants, and trapped-surface diagnostics.

## Independent cross-check and second hierarchy

The 3-D RHS reduces to the separately implemented 1-D ADM RHS when fields
depend on only one coordinate, to roundoff on a 12-point manufactured
state. A flat-space manufactured momentum test and a conformal Ricci
identity test also pass. A second predeclared 16/32/64 hierarchy using
the **same ring seed and run time** has final Hamiltonian RMS residuals
2.31e-5, 9.40e-6, and 2.74e-6, while momentum RMS residuals are
3.58e-6, 3.70e-6, and 1.78e-6. Physical-state differences are
0.01631 and 0.00906, giving order 0.847. The first momentum refinement
still rises, so **the second hierarchy also fails the gate**. The full
numbers are in `reproduce/coupled_dee_ring_3d_secondary_checkpoint.json`.
The higher-resolution trend is promising but does not establish
convergence. Resolving this trial calls for more efficient numerical
geometry and sufficiently resolved grids; neither the original failed
hierarchy nor this second result should be discarded.
