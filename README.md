# Dynamic Equilibrium Events — Toroidal–Helical Transport

[![DOI](https://zenodo.org/badge/1381814823.svg)](https://doi.org/10.5281/zenodo.22998528)

Computational research repository accompanying the theoretical physics work of Gerald Dahlberg Jr.

## Overview

This repository explores computational implementations of Dynamic Equilibrium Events (DEE) in coupled toroidal–helical transport geometries.

The work investigates the evolution of coupled transport structures, including axial propagation, transverse compression, geometric coupling, and dynamic equilibrium behavior.

## Research Areas

- Dynamic Equilibrium Events (DEE)
- Coupled toroidal–helical transport
- Axial propagation
- Transverse compression
- Computational physics
- Python / NumPy
- OpenUSD
- NVIDIA Omniverse
- Scientific visualization and digital-twin workflows

## Computational Direction

The repository is intended to develop reproducible computational demonstrations of the mathematical framework, including:

1. Definition of the toroidal–helical geometry
2. Parameterized transport fields
3. Time-dependent DEE evolution
4. Axial propagation
5. Transverse response and compression
6. Numerical diagnostics
7. OpenUSD scene representation
8. NVIDIA Omniverse visualization

## Associated Research

This computational repository accompanies the following open-access theoretical work:

**Helical Perturbations in Finite-Matter Toroidal–Helical Transport: Nonlinear Distortion, Stability, and Dynamical Recovery**  
Gerald Ted Dahlberg Jr. (2026) — Zenodo preprint  
Publication DOI (v1): [10.5281/zenodo.22074096](https://doi.org/10.5281/zenodo.22074096)  
All-versions DOI: [10.5281/zenodo.22074095](https://doi.org/10.5281/zenodo.22074095)

The paper develops a nonlinear perturbation framework for finite-matter toroidal–helical transport using four explicit control variables: perturbation amplitude \(A\), axial propagation parameter \(k_z\), azimuthal mode number \(m\), and characteristic temporal frequency \(\omega\).

### Research Record

Additional underlying theoretical work is archived through Zenodo with persistent DOI records.

Author: **Gerald Dahlberg Jr.**

ORCID: **0009-0001-0672-7636**

DEE v2.1 DOI: [10.5281/zenodo.21477633](https://doi.org/10.5281/zenodo.21477633)

### Repository archive

Zenodo also archives released versions of this GitHub repository. The DOI badge at the top of this README resolves to the latest archived repository release; manuscript DOIs above identify the associated research publication separately.

## Status

This repository is under active development.

The numerical models and visualizations presented here are computational investigations of the proposed theoretical framework. They should not be interpreted as experimental confirmation or independent physical validation of the theory.

## Technical Scope & Validation Record

The repository includes an explicitly scoped stationary, cylindrically symmetric, linearized weak-field radial module. Its core numerical integration architecture has been checked analytically against the differential equation that the module states it solves, and the automated test suite includes a smooth-source grid-convergence test.

The repository also contains progressively stronger nonlinear validation checkpoints: action-derived DEE matter diagnostics, nonlinear ADM constraint diagnostics, conformally flat constraint-satisfying initial data, and a reduced dynamical initial-data checkpoint with nonzero matter momentum and nonzero extrinsic curvature. These results are explicitly scoped to their stated reduced ansatz and do not establish global nonlinear stability.

### Reduced Fixed-Background Nonlinear DEE Matter Evolution

A subsequent reproducibility checkpoint evolves the nonlinear DEE radial matter field from a finite perturbation about the finite branch on a **fixed local background**. For the recorded test (perturbation amplitude 0.02 through `t=8`), the finite branch remains positive and bounded and the resolution differences decrease consistently with the second-order numerical scheme.

This checkpoint is intentionally limited: **the spacetime metric and extrinsic curvature are not evolved.** It therefore does not establish stability of the coupled Einstein–DEE system, preservation of the Einstein constraints during evolution, global toroidal stability, or survival under gravitational backreaction. The machine-readable record is `reproduce/reduced_dee_evolution_checkpoint.json`, and the executable checkpoint is `reproduce/evolve_dee_reduced.py`.

### First Restricted Coupled Einstein–DEE Evolution

The [3+1 source derivation](docs/DEE_3PLUS1_SOURCE_PROJECTION.md) implements and independently checks the full nonlinear matter projections, including the mixed Theta–Psi kinetic term. A [short plane-symmetric ADM evolution](docs/COUPLED_ADM_1D_CHECKPOINT.md) then evolves the diagonal spatial metric, extrinsic curvature, and all three DEE scalars together on a periodic domain. Across 64, 128, and 256 grid points through `t=0.12`, the Hamiltonian and momentum constraint residuals decrease on refinement, with observed state convergence order 1.998. The [executable](reproduce/evolve_dee_adm_1d.py) and [machine-readable results](reproduce/coupled_dee_adm_1d_checkpoint.json) record the calculation.

This is a plane-symmetric, short-duration test with unit lapse and zero shift. It is **not** a toroidal/helical BSSN evolution and makes no claim of a wormhole solution or nonlinear stability. The remaining validation gate is a full geometric formulation and gauge appropriate to the target geometry, constraint-satisfying data in that geometry, and independently monitored longer three-resolution evolution including curvature and trapped-surface diagnostics.

A further [symmetry-restricted BSSN checkpoint](docs/COUPLED_BSSN_1D_CHECKPOINT.md) evolves the conformal metric, traceless curvature, trace, conformal factor, and contracted connection with the same nonlinear DEE source. Its conformal Ricci operator uses the evolved connection. At 64, 128, and 256 points, the Hamiltonian, momentum, and connection constraints decrease with refinement; the physical fields also approach the independent ADM implementation at approximately second order. This remains a one-dimensional diagonal test with fixed gauge, rather than the target multidimensional toroidal/helical evolution. The [executable](reproduce/evolve_dee_bssn_variables_1d.py) and [numerical record](reproduce/coupled_dee_bssn_variables_1d_checkpoint.json) are provided for reproduction.

An adversarial AI-assisted technical audit conducted on 2026-09-24 examined the radial BVP formulation, finite-domain normalization, repository scope, stress-energy assumptions, and candidate tensor-aware treatment of angular/twist sectors. During that audit, earlier claims that the repository had been demonstrated to violate stress-energy conservation and that it exhibited a "fatal mathematical failure" were withdrawn after comparison with the actual implementation and its documented scope.

The audit does **not** constitute independent peer review or experimental validation.

The complete technical record is archived here:

**[Adversarial AI Technical Audit — DEE Weak-Field Module](docs/AI_TECHNICAL_AUDIT_2026-09-24.md)**

## License

Software in this repository is released under the MIT License unless otherwise noted.

## Reproducible computational companion

The repository contains executable numerical checkpoints and automated tests. Install the dependencies declared by the repository and run the test suite from the repository root.

### Minimal finite-equilibrium branch

`dee_field_model.py` implements the algebraic finite-equilibrium branch

\[
\rho_*^2 = v^2 - \frac{aX_\Theta+bX_\Psi+2cX_{\Theta\Psi}}{\lambda_V},
\qquad m_{\rm eff}^2=2\lambda_V\rho_*^2.
\]

At a supplied field value \(\rho^2\), it uses \(K_\Theta=1+a\rho^2\), \(K_\Psi=1+b\rho^2\), and \(\Lambda=c\rho^2\).

### NEC diagnostic

`nec_scan.py` evaluates the supplied quadratic-form diagnostic \(T_{kk}=u^T K u\). A finite numerical scan is distinguished from the corresponding analytic positive-definiteness condition.

### Linearized cylindrical radial calculation

`linearized_einstein.py` numerically integrates the explicitly scoped stationary cylindrical linearized weak-field radial equation. It is not represented as a nonlinear general-relativity solver.

## Scope and limitations

The included relations and computations are theoretical and numerical tools for exploring proposed DEE assumptions. They do **not** provide experimental validation, a proof of global non-collapse, or a completed full nonlinear Einstein–DEE evolution. Boundary data, source profiles, parameter choices, gauge/ansatz choices, and the interpretation of transport variables remain modeling inputs and hypothetical aspects of the framework.
