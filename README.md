# Dynamic Equilibrium Events — Toroidal–Helical Transport

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
DOI: [10.5281/zenodo.22074096](https://doi.org/10.5281/zenodo.22074096)

The paper develops a nonlinear perturbation framework for finite-matter toroidal–helical transport using four explicit control variables: perturbation amplitude \(A\), axial propagation parameter \(k_z\), azimuthal mode number \(m\), and characteristic temporal frequency \(\omega\).

### Research Record

Additional underlying theoretical work is archived through Zenodo with persistent DOI records.

Author: **Gerald Dahlberg Jr.**

ORCID: **0009-0001-0672-7636**

DEE v2.1 DOI: [10.5281/zenodo.21477633](https://doi.org/10.5281/zenodo.21477633)

## Status

This repository is under active development.

The numerical models and visualizations presented here are computational investigations of the proposed theoretical framework. They should not be interpreted as experimental confirmation or independent physical validation of the theory.

## License

Software in this repository is released under the MIT License unless otherwise noted.

## Reproducible computational companion

The scripts in this repository use only [NumPy](https://numpy.org/) and the
Python standard library. Install the sole numerical dependency and run the automated checks from the repository root:

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

### Minimal finite-equilibrium branch

`dee_field_model.py` implements the algebraic finite-equilibrium branch

\[
\rho_*^2 = v^2 - \frac{aX_\Theta+bX_\Psi+2cX_{\Theta\Psi}}{\lambda_V},
\qquad m_{\rm eff}^2=2\lambda_V\rho_*^2.
\]

At a supplied field value \(\rho^2\), it uses
\(K_\Theta=1+a\rho^2\), \(K_\Psi=1+b\rho^2\), and
\(\Lambda=c\rho^2\), giving
\(K=\begin{pmatrix}K_\Theta&\Lambda\\\Lambda&K_\Psi\end{pmatrix}\). It computes
\(\det K=1+(a+b)\rho^2+(ab-c^2)\rho^4\), and applies the
2-by-2 Sylvester criterion \(K_\Theta>0\) and \(\det K>0\) for
positive definiteness. These are direct
algebraic consequences of the definitions encoded in the program; they do not
by themselves establish the physical applicability of a DEE model.

### NEC diagnostic

`nec_scan.py` evaluates the supplied quadratic-form diagnostic
\(T_{kk}=u^T K u\). Its numerical scan samples unit two-component transport
projections \(u=(\cos\phi,\sin\phi)\). A scan is a finite numerical
sampling diagnostic; positive definiteness of the specified 2-by-2 matrix is
the corresponding analytic all-directions condition within this minimal model.

### Linearized cylindrical radial calculation

`linearized_einstein.py` numerically integrates the stationary radial equation

\[
\frac{1}{r}\frac{d}{dr}\left(r\frac{d\bar h_{\mu\nu}}{dr}\right)
= -16\pi G T_{\mu\nu}(r).
\]

The integration enforces regular-axis behavior through
\(r\,d\bar h/dr\to0\) at \(r=0\), and it accepts a finite outer value for
\(\bar h_{\mu\nu}\). It can solve one radial source profile or multiple
trailing tensor-component profiles independently. Each call emits a warning:
this is a **stationary, cylindrical, linearized weak-field approximation**, not
a nonlinear general-relativity solution.

### Tests

`tests/test_dee_companion.py` verifies the analytic equilibrium and effective
mass expressions, kinetic-matrix positivity, NEC quadratic forms and scans,
and second-order grid convergence of the radial integration for a smooth source.

## Scope and limitations

The included relations and computations are theoretical and numerical tools for
exploring proposed DEE assumptions. They do **not** provide experimental
validation, nonlinear stability, a proof of global non-collapse, or a full
nonlinear Einstein-equation solution. Boundary data, source profiles, parameter
choices, and the interpretation of transport variables remain modeling inputs
and hypothetical aspects of the framework.
