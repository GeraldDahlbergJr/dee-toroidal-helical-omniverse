# Single manufactured qualification protocol

Base 5a83752a0e8df98e9db5c25041aebcd3e87ed6f9; branch
run43-independent-physical-validation. Wall budget starts 2026-11-06
14:57:04 UTC, deadline 15:27:04 UTC. One logical CPU, BLAS/OpenMP one
thread, 4 GiB address-space limit. Only n=8,16,32. No physical-source solve.

## Representation and weak equations

There are n uniform spans in each coordinate. Degree 4, simple interior
knots, open radial end knots of multiplicity 5; remove first/last radial
basis functions. Angular cardinal splines on extended uniform knots are
summed by indices modulo n. Scalar unknown is psi=1+c; three W components
are Cartesian. Thus c and W have identically zero radial traces; angular
periodicity is an identity of functions, including derivatives through order
3. Scalar dimension (n+2)n^2; total coupled dimension 4(n+2)n^2.

Exact chart X=((1+r cos t)cos p,(1+r cos t)sin p,r sin t),
r in [.05,.25], q=1+r cos t, J=rq. Chart reciprocal frame is
E=(er,et/r,ep/q); gradient components are sum_a E_ai d_a.
Volume integration uses J times parameter Gauss weights, no geometry
interpolation. Primary assembly: 6 Gauss points per span per direction.
This rule is not claimed exact for geometric rational/trigonometric factors.

Use the original source A=.0127638,
B=.068614375(r^-2+9/q^2), S=-.030316(et/r-3ep/q), G=.001.
Define N(u,W)=u^-7 |LW|^2/8+2 pi G(Au^5+Bu).
Manufactured equations are -Delta u=N(u,W)+fH,
-div LW=-8 pi G u^6 S+fM. Weak equations are
(grad(u-1),grad v)=(N+fH,v) and (LW,Lz)/2=(fM-8 pi G u^6 S,z).
Assembly explicitly forms the full Cartesian STF strain, including both
symmetric off-diagonal entries. Homogeneous test traces remove boundary
terms; no traction condition is imposed. No stabilization is present.

CG uses tensor generalized-eigenfunction Poisson preconditioning only; this
does not alter the equations. Scalar Newton steps and alternating vector
updates solve the nonlinear fixture. Each linear solve explicitly checks
Euclidean coefficient residual norm ||Ax-b||2/||b||2 <=1e-10
(requested CG tolerance 2e-12). These are algebraic checks, not validation.

## One nonpolynomial analytic fixture

Let f=sin(pi(r-.05)/.20) exp(sin(t)/5+cos(p)/7).

    u=1+.01 f
    Wx=.001 f (1+sin(t+p))
    Wy=.0007 f (1+cos(t-2p))
    Wz=.0009 f (1+sin(2t-p)).

These fields are smooth, periodic, u positive and have prescribed radial
traces. Define Cartesian D_i=E_ai partial_a and
Delta v=v_rr+(1/r+cos(t)/q)v_r+v_tt/r^2-sin(t)v_t/(rq)+v_pp/q^2.
The independently derived forcing is

    fH=-Delta u-u^-7 sum_ij(Dj Wi+Di Wj-2 delta_ij divW/3)^2/8
       -2 pi G(Au^5+Bu)
    fMi=-Delta Wi-Di(sum_j Dj Wj)/3+8 pi G u^6 Si.

SymPy differentiates these analytic expressions directly before numerical
assembly; no spline matrix or discrete derivative enters forcing construction.
The script contains the full expressions as an executable symbolic definition.

## Independent residuals and scoring

Saved coefficients define the approximate continuous fields. A separate module
converts each scalar spline to local power polynomials and evaluates derivatives
by differentiating coefficients and Horner evaluation. It uses independently
coded Cartesian second derivatives

 D_j D_i v = sum_ab E_aj E_bi v_ab + E_aj (partial_a E_bi) v_b.

Then div LW=Delta W+grad(div W)/3. It does not use the weak operator,
Space.grad, or a composed discrete first derivative. Physical-scaled fixture
residuals are H=-8 u^-5(Delta u+N+fH),
Mi=u^-10(div LW_i-8 pi G u^6 Si+fMi).
Report normalized volume-weighted RMS with J u^6 in inner (.05,.10),
bulk [.10,.20], outer (.20,.25), plus surfaces r=.05+.20/n and
r=.25-.20/n (first interior planes); surfaces use J u^6 angular weights.
The 6- and 8-point per-span integration rules are independent from the solve
when applied to the power-polynomial evaluator. RMS field error uses the
same volume weight; maximum field error is sampled, not a supremum bound.

Every component and region must have order >=1.5 on both pairs. Quadrature
uncertainty proxy uE=2|E8-E6| gives conservative-for-that-proxy
p_low=log2((Ec-uc)/(Ef+uf)). It is not a certified bound on other error
sources. Total evaluator, quadrature, and solver uncertainty must be assessed
before qualification; diagnostic agreement is not silently promoted to a bound.

Trace norm: maximum absolute field defect, with declared characteristic scales
[1,.0025,.00175,.00225] for [u,Wx,Wy,Wz]. Thresholds are 100 eps times
these fixed scales (eps is IEEE float64). Strong construction gives an entire
trace identity; samples independently check implementation roundoff. No scale
is changed in response to results. Compare derivative evaluation against the
assembly spline implementation through order 2. No unexpected discrete kernel
is permitted; examine basis independence, coercivity and local coupled effects.

Stop immediately on a decisive criterion failure or resource exhaustion.
Missing evidence is reported as unresolved and cannot produce QUALIFIED.
No further correction candidate or original-source solve is authorized.
