# Stage A½: recovered inputs and exact symbolic verification

Date: 2026-11-10 UTC. Base: `bd32e856b49bc0c227d235ef68f345d72acb2b21`.
Repository: GeraldDahlbergJr/dee-toroidal-helical-omniverse.
Branch: `research/stage-a-half-symbolic-recovery`.

## Source recovery and scope

The analytical candidates were supplied explicitly by the author in the current research request. They are newly available inputs, not recovered historical files, an authenticated transcript, or evidence for the missing historical Gaussian source/generator/configuration. No missing historical commit is reused or reconstructed.

The action is present in `evidence/einstein_dee_sources/manuscripts/DEE_LaTeX_v2/main.tex`, Eqs. (3)–(9), (14)–(15), and transcribed in `docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md`. Put Q=(Theta,Psi)^T and A=[[a,c],[c,b]]. Its matter part is

    L_m = -1/2 g^{mu nu} partial_mu rho partial_nu rho - V(rho)
          -1/2 g^{mu nu} (partial_mu Q)^T M(rho) partial_nu Q,
    M(rho)=I+rho^2 A,  V(rho)=lambda_V (rho^2-v^2)^2/4.

Use Cartesian (t,x,y,z), eta=diag(-1,1,1,1), Box=-partial_t²+partial_x²+partial_y²+partial_z². All constants are spacetime independent. lambda_V>0, v>0, rho_0>0 and rho_0²=v²-C/lambda_V>0. A is real symmetric; k is a real two-vector. Positive-definite M(rho_0) is an additional physical admissibility assumption inherited from the matter-existence result, not necessary for these algebraic identities. No parameter values are fitted. Exact polynomial manipulations substitute v²=rho_0²+C/lambda_V; physical v²>0 must also hold when selecting parameters.

Define s=x²+y², C=k^T A k, B=k^T M(rho_0)k, epsilon=V(rho_0)+B/2. The supplied configuration rho=rho_0, Q=kz+Q_0 gives

    E_rho^(0)=-rho_0[lambda_V(rho_0²-v²)+C]=0,
    E_Q^(0)=(0,0)^T,
    T^(0)=diag(epsilon,-epsilon,-epsilon,B-epsilon).

Both phase currents and this stress are constant, with zero flat divergence. This recovers an extended matter configuration, not a localized source.

## Metric calculation

The supplied trace-reversed perturbation has the following explicit nonzero components:

    bar h_tt = -4 pi G epsilon s
    bar h_zz = -4 pi G (B-epsilon) s
    bar h_xx = 6 pi G epsilon s - 4 pi G epsilon x²
    bar h_yy = 6 pi G epsilon s - 4 pi G epsilon y²
    bar h_xy = bar h_yx = -4 pi G epsilon xy.

Each residual R_mu nu=Box bar h_mu nu+16 pi G T_mu nu^(0) is checked separately:

| Component | Box bar h | Source -16 pi G T^(0) | Residual |
|---|---|---|---|
| tt | -16 pi G epsilon | -16 pi G epsilon | 0 |
| tx | 0 | 0 | 0 |
| ty | 0 | 0 | 0 |
| tz | 0 | 0 | 0 |
| xx | 16 pi G epsilon | 16 pi G epsilon | 0 |
| xy | 0 | 0 | 0 |
| xz | 0 | 0 | 0 |
| yy | 16 pi G epsilon | 16 pi G epsilon | 0 |
| yz | 0 | 0 | 0 |
| zz | -16 pi G (B-epsilon) | -16 pi G (B-epsilon) | 0 |

Gauge residuals partial^mu bar h_mu nu are (0,0,0,0). Explicitly the x constraint is 4 pi G epsilon x-4 pi G epsilon x=0 and the y constraint is -4 pi G epsilon y+4 pi G epsilon y=0; t,z constraints vanish by time and axial independence. The transverse Cartesian tensor terms are essential to these cancellations.

Independent trace reversal gives

    bar h = eta^{mu nu} bar h_mu nu = 4 pi G (4epsilon-B)s,
    h_mu nu = bar h_mu nu - (1/2) eta_mu nu bar h,
    h_tt = 2 pi G (2epsilon-B)s,
    h_xx = 2 pi G (B-epsilon)s - 4 pi G epsilon x²,
    h_yy = 2 pi G (B-epsilon)s - 4 pi G epsilon y²,
    h_xy = -4 pi G epsilon xy,
    h_zz = -2 pi G (B+2epsilon)s,
    h = 4 pi G (B-4epsilon)s.

All other components vanish. The verifier also constructs linearized Ricci directly from h, without using the wave-equation shortcut, and obtains all sixteen entries of G^(1)-8 pi G T^(0) equal to zero. Trace reversal applied twice reproduces bar h exactly.

## Claim classification

| Claim | Classification | Qualification |
|---|---|---|
| Supplied constant matter configuration solves the flat equations and gives the prescribed stress | VERIFIED | Branch relation and original action as above |
| Candidate solves all ten sourced metric equations | VERIFIED | Local linearized calculation |
| Candidate satisfies all four harmonic constraints | VERIFIED | Cartesian tensor calculation |
| Ordinary h and h_zz above | VERIFIED | Trace reversal and direct Einstein-tensor check |
| General amplitude and fixed-phase transverse equations | VERIFIED | Derivation and signs in companion backreaction check |
| Reported frozen amplitude residual | VERIFIED | Defined as the full amplitude EOM left side, or reduced LHS minus RHS |
| Frozen matter supplies a generic first-order matter solution with this h | CORRECTION REQUIRED | Nonzero amplitude residual unless its coefficient vanishes |
| First-order conservation proves the amplitude equation | CORRECTION REQUIRED | Constant background rho makes that implication invalid |
| This candidate is a global weak-field, localized or complete coupled spacetime | NOT YET ESTABLISHED | No global boundary/matching solution or uniform weak-field bound |
| Historical provenance beyond the explicitly supplied inputs | NOT YET ESTABLISHED | Missing historical artifacts remain missing |

These are exact algebra checks by the same investigator using distinct calculation routes, not external peer review. The perturbation grows as s. Where G times the stress scales times s is not small, a linear truncation cannot establish validity. No claim about global asymptotics, nonlinear stability, regular complete spacetime, or localized gravitational matching follows.

## Reproduce and preserve

From the repository root, with Python 3 and SymPy (checked with 1.14.0):

```sh
python DEE_STAGE_A_HALF_VERIFIER.py
```

The verifier prints each of the ten equation residuals and four gauge residuals as zero, the ordinary h, action-derived matter/stress checks, four conservation identities, frozen residuals and a final PASS. It raises an exception for any failed equality. The exact frozen amplitude residual is nonzero in general; PASS checks its value, not that frozen matter solves the equation.

Standing verdicts from `DEE_SCIENTIFIC_STATUS.md` remain unchanged: **MATTER EXISTENCE: ESTABLISHED** for the previously delimited extended flat-background family; **LOCALIZATION: OBSTRUCTION IDENTIFIED** under the existing theorem's stated hypotheses; **COMPLETE EINSTEIN–DEE SPACETIME: NOT YET ESTABLISHED**. The formal Eq. (34) correction remains required. No pre-existing evidence is edited. No Run 43 rerun, numerical PDE solution, parameter fitting, speculative model change, or nonlinear gravitational construction is performed.
