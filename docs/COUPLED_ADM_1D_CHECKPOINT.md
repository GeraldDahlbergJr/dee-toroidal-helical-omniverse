# First restricted coupled Einstein–DEE evolution

Run `PYTHONPATH=. python reproduce/evolve_dee_adm_1d.py` from the repository root.
The reproducible numerical output is in
`reproduce/coupled_dee_adm_1d_checkpoint.json`.

This is a plane-symmetric periodic 3+1 **ADM** evolution with diagonal
`gamma_ij` and `K_ij`, unit lapse, zero shift, and all three nonlinear
action-derived DEE fields. It uses centered second-order spatial differences
and RK4 time stepping to `t=0.12`. The transport fields start as spatial
sinusoids with nonzero mixed gradient energy (`c=0.1`); their momenta and the
radial scalar respond during evolution. The periodic Lichnerowicz equation
supplies approximately constraint-satisfying initial data with a spatially
constant isotropic mean curvature; the initial residual is due to evaluating
the spectral initial solve with the evolution's finite-difference operator.
The source is the independently checked `E`, `S_i`, and `S_ij` projection.

| Points | Initial H RMS | Final H RMS | Final M_x RMS |
| ---: | ---: | ---: | ---: |
| 64 | 8.8782e-7 | 8.5796e-7 | 1.1334e-7 |
| 128 | 2.2335e-7 | 2.1580e-7 | 2.8596e-8 |
| 256 | 5.5925e-8 | 5.4031e-8 | 7.1653e-9 |

Adjacent-grid maximum state differences are `6.2024e-6` and `1.5525e-6`,
giving an observed order of `1.9983`. All evolved spatial metric components
stay positive, the transport kinetic determinant stays positive, and
`rho_min > 0` over this short interval. The diagnostic evaluates the
Hamiltonian `R+K^2-K_ij K^ij-16*pi*G*E` and momentum
`D_j K^j_i-D_i K-8*pi*G*S_i` independently from evolved fields.

The symmetry restriction excludes toroidal/helical angular geometry,
off-diagonal metric components, gravitational radiation, boundaries, and
trapped-surface inference. Geodesic slicing has no long-term stability
guarantee. This run demonstrates a short coupled local evolution and
convergent constraints for these data; it does not establish a BSSN
implementation or stability of a wormhole or full DEE configuration.
