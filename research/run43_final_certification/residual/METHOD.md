# Residual workstream: fixed-field audit, not certification

The coefficient archives and every entry in the completion SHA256SUMS were
checked. `reproduce.py` evaluates archived fields, never calls `solve`, and
writes all 60 scores for rule selected by CERT_RULE (default 8). Binary64
agreement is reproduction, not an enclosure. All 60 q6 and all 60 q8 scores
reproduced bit identically. Cumulative q6 evaluation times at n=8/16/32 were
0.476/2.725/17.427 s; q8 times were 0.798/5.090/35.963 s. Every manifest
hash matched. Outputs are reproduction_q6.json and reproduction_q8.json.

## Pilot calculations

`pilot.py` uses exact Fraction arithmetic on the archived binary64 score.
For n=32 first_outer My, E=5.722436413449492e-7 and its total allowance
is 5.722436413449492e-8. The q6/q8 difference is
3.6089171206474615e-14. This difference does not bound truncation.

For a C^16 function on each angular cell, the Gauss-8 remainder constant
for an interval of length h is C8*h^17 sup|f^(16)|, where
C8=(8!)^4/[17(16!)^3]=1/58914962958794342400000.
Telescoping the tensor product with positive weights gives surface error
at most L^2 C8 h^16 (sup|d_t^16 F|+sup|d_p^16 F|), where
F=r*q*u^6*My^2, L=2*pi, h=L/32. Piecewise derivatives suffice because
quadrature cells follow spline cells. At r=39/160, u>=1, q>=3/4,
the denominator is >=1053/160 using pi>=3. Allocating only
0.10*E^2*Dlower=2.1551144541743275e-13 to numerator truncation would
be sufficient if the derivative sum were <=6.542024184414019e19
(using pi<=22/7). This is a **target**, not an established derivative bound;
denominator truncation, arithmetic, and solution error require additional
allocations. No numerical derivative samples can certify that supremum.

## A valid but unusably coarse continuum enclosure

`coarse_pilot.py` supplies a pointwise analytic bound, bypassing quadrature.
Interpret the stored binary64 coefficients as exact reals. The script checks
|c|<=1 and c_scalar>=0; positivity and partition of unity give 1<=u<=2.
It checks every positive binary64 knot gap is >=h0=1/(10*n). The exact
spline derivative recurrence bounds each parameter first derivative by
D1=8/h0 and each pure/mixed second derivative by D2=64/h0^2. Missing
radial boundary basis coefficients are zero, so the same bounds apply.
Across r in [1/20,1/4], q>=3/4, each frame entry |E_ai|<=22 and each
frame coordinate derivative |d_a E_bi|<=500. Consequently each Cartesian
Hessian entry is bounded by
H=9*(22^2 D2+22*500 D1)=4356 D2+99000 D1.
The vector operator is lap(w_i)+(1/3) sum_j d_i d_j w_j,
whose component magnitude is <=4H.

For the exact manufactured fields, all parameter first derivatives are
<=100 and all second derivatives <=10000. These deliberately loose
bounds follow by the product rule: exp((1/5)sin(t)+(1/7)cos(p))<2;
radial frequency 5*pi<16; angular exponent frequencies <=1; modulation
frequencies <=2; field amplitudes <=1/100. Even the loose common bound
1 for each field amplitude gives first derivative <=2*16+4<100,
and second derivative <=2*16^2+2*16*4+12<10000.
Thus each exact Cartesian Hessian entry is <=Hex=4356*10000+99000*100,
and each exact vector operator component <=4Hex. Both exact and stored
scalar fields are >=1 and <=2. The source |8*pi*.001*u^6*S_i| is below
200000 (indeed below 3), since |S_i|<.031*(20+4)<1. Therefore
|My|<=4H+4Hex+400000=115417926400. Any positive weighted RMS over
the pilot surface belongs to [0,115417926400]. Fraction arithmetic gives
an absolute difference bound from E of 115417926400 (rounded display).
This is a rigorous fixed-field range bound, over 2e18 times the allowance,
and explicitly **does not** imply failure of the actual score or include
an exact-assembly solution neighborhood. The useful pilot test failed to
establish an adequately sharp bound.

The binary knot range and ideal [1/20,1/4] differ by tiny endpoint rounding.
For the pilot r=39/160 lies strictly inside both. All angular arguments and
radii used in this proof are continuum real variables, not computed nodes;
the intentionally loose frame constants cover the binary endpoint range too.

## Precise remaining obstruction

The continuum bound destroys the cancellation between saved-field second
derivatives and manufactured source. Sharp certification needs local,
cancellation-preserving bounds for F and its derivatives, or a rigorous
interval integral enclosure. The existing PPoly coefficient conversion,
quadrature nodes/weights, transcendental evaluation, and reductions have
no verified reference enclosure. No bound for these operations follows
from replay matching the archive. Surface integration alone also cannot
supply propagated solution error. The 16th derivative ceiling above is a
concrete next test for a cellwise Taylor/interval implementation; it must
be established on every angular cell, together with denominator and
arithmetic budgets, before extending to 60 scores.

Reproduce from repository root:

    CERT_RULE=8 OPENBLAS_NUM_THREADS=1 python research/run43_final_certification/residual/reproduce.py
    CERT_RULE=6 OPENBLAS_NUM_THREADS=1 python research/run43_final_certification/residual/reproduce.py
    python research/run43_final_certification/residual/pilot.py
    python research/run43_final_certification/residual/coarse_pilot.py

Verdict remains NOT YET ESTABLISHED.
