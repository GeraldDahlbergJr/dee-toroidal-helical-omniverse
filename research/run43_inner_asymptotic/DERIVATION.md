# Consistency of the frozen toroidal momentum operator

## Operator and radial rows

Use Cartesian components throughout. At fixed angles let $D$ be the frozen order-two radial first derivative, with uniform spacing $h$. Its inner wall stencil is $(-3,4,-1)/(2h)$ on offsets $0,1,2$. Interior rows use $(-1,0,1)/(2h)$. The toroidal gradient uses $A_r=e_r$, $A_t=e_\theta/r$, $A_p=e_\phi/q$, where $q=R+r\cos\theta$. For vectors $v,a$, write
$$
L(v,a)=v\otimes a+a\otimes v-\tfrac23(v\cdot a)I.
$$
The momentum divergence is the sum, over outer direction $i$ and inner direction $j$, of $[\partial_i L(\partial_j W,A_j)] A_i$. The frozen implementation substitutes $D_i,D_j$, using periodic centered first differences in both angles. Both evaluate $\psi^{-10}$ times the divergence minus $8\pi G\psi^{-4}S$.

The following Taylor expansions follow from the stencil moments about **each row's own node**. The zeroth and first moments vanish and the second moment of the dimensionless composed stencil equals 2.

| Radial row | Offsets | Coefficients of $h^2D^2$ | $D^2W-W_{rr}$ |
| --- | --- | --- | --- |
| Inner wall $0$ | $0,1,2,3$ | $(5,-11,7,-1)/4$ | $3hW_{rrr}/4+5h^2W_{rrrr}/24-h^3 W_{rrrrr}/16+O(h^4)$ |
| First interior $1$ | $-1,0,1,2$ | $(3,-5,1,1)/4$ | $hW_{rrr}/4+5h^2W_{rrrr}/24+h^3 W_{rrrrr}/16+O(h^4)$ |
| $2\le k\le N-2$ | $-2,0,2$ | $(1,-2,1)/4$ | $h^2W_{rrrr}/3+2h^4W_{rrrrrr}/45+O(h^6)$ |

The outer wall and its neighbor reflect the inner stencils, reversing odd terms. For example, the first interior composition differentiates the centered first derivative at row 2 and the one-sided derivative at row 0. Their leading first-derivative error coefficients differ by $1/6-(-1/3)=1/2$, leaving $hW_{rrr}/4$. The wall results require one-sided smooth traces through the stencil support. The displayed higher remainders use bounded derivatives through order 6 and 8 as applicable. The leading-order conclusions require only enough uniformly bounded radial and angular derivatives to control the stated lower-order remainders.

Since $A_r$ is independent of $r$, the pure radial term has no metric commutator. Its error in momentum is $\psi^{-10} B_r(D^2W-W_{rr})$, with $B_r=I+e_r e_r^T/3$. This map has eigenvalues $1,1,4/3$, so it preserves nonzero third radial traces in the defect.

## Metric terms and cancellation

For a scalar coordinate direction $i$, put $v_h=D_jW$ and $A=A_j$. The exact discrete product-rule remainder before contraction with $A_i$ is
$$
 C_i(v_h,A)=D_i L(v_h,A)-L(D_i v_h,A)
 =\sum_s d_{i,s}\,L(v_h(x+s\delta_i),A(x+s\delta_i)-A(x)).
$$
Here $d_{i,s}$ are the first-derivative coefficients including $1/\delta_i$. For a smooth $v=\partial_jW$ and a centered outer stencil,
$$
 C_i(v,A)-L(v,\partial_i A)=
 \tfrac{\delta_i^2}{6}
 [L(v,\partial_i^3 A)+3L(\partial_i v,\partial_i^2 A)
       +3L(\partial_i^2 v,\partial_i A)]+O(\delta_i^4).
$$
At the radial wall the same third-derivative combination has coefficient $-h^2/3$ with $O(h^3)$ remainder. In both formulas one must add $C_i(v_h-v,A)$. The stencil representation gives a bound $O(\|v_h-v\|_\infty)=O(\delta_j^2)$ because metric differences are $O(\delta_i)$. On smooth stencils the leading contribution from this part is the inner first-derivative error multiplied by $\partial_i A$. This accounts for both error sources in a metric defect measured against a continuum derivative.

In particular $\partial_r A_t=-e_\theta/r^2$, $\partial_r^2 A_t=2e_\theta/r^3$, $\partial_r^3 A_t=-6e_\theta/r^4$, and $\partial_\theta A_r=e_\theta$, $\partial_\theta^2 A_r=-e_r$, $\partial_\theta^3 A_r=-e_\theta$. They enter $rt$ and $tr$, respectively, with different vector contractions. There is no general identity canceling their errors. Terms with opposite signs can cancel in the vector sum, including $tt$ and the remaining principal terms. The mixed principal terms are second order for smooth data, since a derivative in a distinct coordinate preserves the order of the inner first difference. Pure angular compositions are second order on the periodic, smooth fields considered here. Consequently no metric or mixed term cancels the order-$h$ radial truncation coefficient asymptotically under simultaneous refinement with angular steps $O(h)$, although finite-grid cancellations can be large.

## Fixed-width weighted RMS

For the gate, $a=0.05$, $b=0.10$, $h=0.20/N$, and the wall is excluded. There are $m=\lfloor(b-a)/h\rfloor$ radial rows, exactly $N/4$ in these experiments. The weight is $\rho=rq\psi^6$, positive and bounded above and below in the band if $\psi$ is bounded and bounded away from zero. The weight contains neither a singular axis nor a vanishing wall factor. Multiplying numerator and denominator by $h\,\Delta\theta\,\Delta\phi$ leaves the frozen RMS unchanged. Uniform angular counts cancel in the first-row fraction even for unequal poloidal and azimuthal counts. For smooth periodic fields the row's angular weighted energy converges by quadrature, and
$$
 \frac{\sum_{\theta,\phi}\rho(a+h)}{\sum_{k=1}^m\sum_{\theta,\phi}\rho(a+kh)}
 = h\frac{\int \rho(a)\,d\theta d\phi}
              {\int_a^b\int \rho\,d\theta d\phi\,dr}+O(h^2)
$$
on the co-refined grids with smooth data. The $O(h^2)$ includes the endpoint quadrature error of the radial rectangle rule. Let $\tau_h$ be the frozen minus continuum operator on exact samples and define $A=\psi(a)^{-10}B_r W_{rrr}(a)/4$. Uniform radial Taylor expansions up to the wall and angular resolution with second-difference error $O(h^2)$ give $\tau_h(a+h)=h A+O(h^2)$, and uniformly $O(h^2)$ on the remaining rows. Thus
$$
 \|\tau_h\|_{\rm band}^2
 = C h^3+O(h^4),\qquad
 C=\frac{\int \rho(a)|A|^2\,d\theta d\phi}
         {\int_a^b\int \rho\,d\theta d\phi\,dr}.
$$
The cross term between $hA$ and a same-row $O(h^2)$ error is $O(h^4)$ in band-normalized energy and may have either sign. Separate radial rows have disjoint support and no cross term. The remaining band contributes nonnegative $O(h^4)$ energy. If $C>0$, the asymptotic order is exactly $3/2$. If $C=0$ under these smoothness hypotheses, the norm is $O(h^2)$, with second order if the bulk coefficient is nonzero. A degenerate wall weight, a band containing only a fixed number of rows, unresolved angular structure, or nonsmooth boundary traces require new estimates. In particular a slower angular refinement can make the angular error compete with or dominate $h^{3/2}$.

With a regular next-order expansion $E_h^2=C h^3(1+d h+o(h))$ and $C>0$, a doubling pair gives $p_h=1.5+d h/(4\ln 2)+o(h)$. The sign of $d$ depends on radial quadrature, trace variation, same-row cancellation and bulk energy. Even when $p_h>1.5$ for individual pairs, its positive margin vanishes as $h\to0$. Without a smooth expansion for the correction, the leading equivalence still yields $p_h\to1.5$.

## Numerical solutions and the threshold

For exact field samples $U^*$, consistency compares $M_h(I_h U^*)$ with $I_h M(U^*)$; this is the truncation defect above. The frozen residual of computed fields $U_h$ is the sum of this defect, the exact continuum residual (zero for an exact solution), and $M_h(U_h)-M_h(I_h U^*)$. Small nodal error, or a small algebraic residual in a different discretization, provides no automatic bound on the last term, which contains differentiated field errors. The $3/2$ result transfers to computed fields if this last term is $o(h^{3/2})$ in band RMS, or if a row/bulk expansion establishes the same nonzero leading coefficient in the **combined** error. A mesh-dependent error can cancel the row defect or add a larger contribution. A solver correction would have to control that differentiated error while retaining independent continuum convergence; the stencil proof alone supplies no such construction.

The frozen operator is consistent for smooth fields on this shell, first order locally on the wall and first interior rows in general, and second order away from those rows on co-refined grids. The positive first-row coefficient establishes a conditional upper bound on the asymptotic order in the specified RMS. A universal convergence barrier for Run 43 numerical fields does not follow.
