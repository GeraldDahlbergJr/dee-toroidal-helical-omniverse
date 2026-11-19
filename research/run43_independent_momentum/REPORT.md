# Run 43 — independent validation of compatible momentum

**Decision B — Partial improvement.** A conforming weak momentum operator on the original shell substantially reduces independently reconstructed physical momentum residuals at fixed archived conformal factors. The first-row and endpoint discrepancies shrink across all four grids, and mixed-term differences under refinement shrink. Pure angular accuracy does not show a clear improvement. Coupled Hamiltonian–momentum solution accuracy and pointwise boundary consistency of the native finite element field remain open.

## Physical problem and derivation

Use the [authoritative Run 43 specification](../RUN43_CONSOLIDATED_MATHEMATICAL_SPECIFICATION_20261106.md). The exact map is
$X=((1+r\cos\theta)\cos\phi,(1+r\cos\theta)\sin\phi,r\sin\theta)$.
The shell is $0.05\le r\le0.25$, with both angles periodic, $q=1+r\cos\theta$, and $dV=rq\,dr\,d\theta\,d\phi$. The toroidal background metric is $dr^2+r^2d\theta^2+q^2d\phi^2$, which is Euclidean in Cartesian components. For the coordinate gradient introduce $b^r=e_r, b^\theta=e_\theta/r, b^\phi=e_\phi/q$. The physical metric is $\gamma=\psi^4\bar\gamma$ and $A_{\rm phys}^{ij}=\psi^{-10}(LW)^{ij}$.

The original seed has $G=.001$, $\rho^2=.945$, $\Pi=(0,.12,-.04)$, and phase gradients $(.25,.20)d(\theta-3\phi)$. Raising the fixed momentum covector with the Euclidean background gives $S=-.030316(b^\theta-3b^\phi)$. The physical energy is $.0127638+.068614375\,\psi^{-4}(r^{-2}+9q^{-2})$. The unchanged Dirichlet traces are $W=0, \psi=1$ at both radial faces.

For $f=8\pi G\psi^6S$, the continuum equation is $\operatorname{div}LW=f$, where $LW=\nabla W+\nabla W^T-\frac23I\operatorname{div}W$. For any symmetric trace-free $T$, integrating in the Euclidean physical volume gives

```text
(div T,z) = -(T,Lz)/2 + <T n,z>_boundary.
```

The radial normals are $-e_r,+e_r$ with face measure $rq\,d\theta\,d\phi$. Periodic seam contributions cancel for smooth periodic tensors, or pair into jumps for broken tensors. The test functions vanish on the prescribed Dirichlet faces. For nonzero test traces the boundary work survives, and the exact strain-form reaction on a lift is $(f,z)+(LW,Lz)/2$. No traction is imposed. The pilot uses the homogeneous-test form only and computes no pointwise traction.

A second integration by parts on homogeneous-trace fields gives $(LW,Lz)/2=(\nabla W,\nabla z)+(\operatorname{div}W,\operatorname{div}z)/3$. This extends by density to periodic $H^1$ fields with zero radial trace. The pilot assembles this last form. On each chart cell the derivatives are taken directly on conforming basis functions. For the full strong equation the Cartesian Hessian coefficient is

```text
H_ij(v) = sum_(a,b) b_i^a [ b_j^b (partial_ab v) + (partial_a b_j^b) (partial_b v) ]
(div LW)_i = sum_j H_jj(W_i) + sum_j H_ij(W_j)/3.
```

This includes radial and both angular second derivatives, derivatives of reciprocal frame coefficients, and both ordered terms in every mixed pair. In the weak form these terms are integrated together through first derivatives and the exact volume geometry. There is no application of a discrete pointwise product rule.

### Implemented discretization and stability

Use continuous tensor-product Q2 in $(r,\theta)$, with two archive nodal intervals per element, periodic poloidal traces, and homogeneous radial elimination. In the toroidal direction express the Cartesian vector as $W=R_z(\phi)u_m(r,\theta)e^{im\phi}$ per Fourier mode. Axisymmetry of the geometry decouples these modes. In the rotated Cartesian frame at $\phi=0$ the rows of the derivative map act in the radial cylindrical, azimuthal, and vertical directions as

```text
cos(theta) d_r - sin(theta)/r d_theta
q^(-1) (i m I + Omega),   Omega(x,y,z)=(-y,x,0)
sin(theta) d_r + cos(theta)/r d_theta
```

Thus even the evaluated $m=0$ mode retains the toroidal derivative of the Cartesian vector. The operator is derived mode by mode for general data. Only $m=0$ was solved, because the original load and conformal-factor archives have that symmetry. The largest discarded deviation was $2.22\times10^{-16}$ in $\psi$ and $1.81\times10^{-19}$ in $W$. The load uses Q2 interpolation of the toroidal mean of the archived $\psi$, evaluated before taking its sixth power. The saved diagnostic pairs use the original $\psi$ arrays byte-for-byte at the array level.

Per unit toroidal angle, let $D$ stack all derivative maps, $C$ denote their divergence contraction, $B$ denote value evaluation, and $H_q$ be positive Gauss weights times $rq$. The negative-divergence stiffness is $K=D^*H_qD+C^*H_qC/3$, with load $-B^*H_q f$. This is a dual operator. Representing it as a weak approximation to $\operatorname{div}L$ in displacement coefficients would give $-M^{-1}K$ with the corresponding mass. There is no nodal identification between this operator and the two finite-difference matrices. There is no stabilization, load modification, interval selection, or coefficient fitting. The exact-integration conforming formulation in the specification uses additional Cartesian polynomial generators and a discontinuous STF space. This pilot uses a simpler conforming primal space and Gaussian integration.

The preconditioner is $P=D^*H_qD$ using a sparse direct factorization, with the original GMRES `rtol=1e-10, atol=1e-13, restart=60, maxiter=100`. The bound $|\operatorname{div}z|^2\le3|\nabla z|^2$ gives $P\preceq K\preceq2P$. Radial and poloidal polynomial derivatives are resolved in the positive quadrature norm. A zero gradient then implies zero displacement by the radial traces. Bounded positive metric factors, exact Gauss integration of squared polynomial chart derivatives and radial Poincaré give mesh-uniform $H^1$ coercivity in a fixed mode. Quadrature consistency of the stiffness follows for fixed mode because all metric factors are smooth on this shell. On each cell their cubic Taylor polynomials times the basis products are integrated exactly by Gauss-4 in each direction, leaving a coefficient remainder bounded by $O(h^4)$ in the derivative and value norms. The nonpolynomial load requires its own consistency estimate. The four-versus-six-point experiment measures sensitivity on this physical input family and supplies no rigorous load-quadrature error enclosure.

## Finite-grid comparisons and independent checks

The original expanded implementation acts as $E_i=L_2W_i+G_i\sum_jG_jW_j/3$, where $L_2$ uses direct chart second differences and $G_i=\sum_a\operatorname{diag}(b_i^a)D_a$. The frozen 12/12 archive instead solves $T W-\frac14[h_\theta^{-2}r^{-2}\delta_\theta^4+h_\phi^{-2}q^{-2}\delta_\phi^4]W$, where
$T_{ij}=\delta_{ij}\sum_kG_kG_k+G_jG_i-2G_iG_j/3$
and $\delta^4$ has coefficients $1,-4,6,-4,1$. Exactly,
$(T-E)_{ij}=\delta_{ij}(\sum_kG_kG_k-L_2)+[G_j,G_i]$.
Both are restricted to interior equations after applying the full stencils with Dirichlet data. The executable applies these exact actions to both sets of samples, compares them to the independent continuum-target reconstruction, and records all interior regional RMS and maxima. The archived fields' verified stabilized relative algebraic residuals are $5.74,9.91,9.88,2.44\times10^{-11}$, respectively. Applying the original expanded operator to those same fields gives relative defects $0.235,0.279,0.230,0.140$. Those are cross-operator defects rather than solve errors of an archived expanded solve.

The first-row analysis gives $(D_rD_r-D_{rr})W_1=(-W_0+3W_1-3W_2+W_3)/(4h_r^2)$, with radial vector factor $I+e_re_r^T/3$ when comparing the tensor action against a direct Hessian. The candidate uses neither second-difference route. The frozen audit still computes composed derivatives on its input arrays. For smooth fields the first-row term is $O(h_r)$, and its coefficient need not be small. A fixed-band volume average can yield asymptotic $O(h_r^{3/2})$ for a single such row, with coarse ratios also affected by near-boundary variation and other terms. The audit mismatch is consistent with the known stencil defect. This run does not attribute its entire frozen residual to that one term. The exact commutator identity retains the angular and mixed coefficient effects, which a radial closure alone cannot address.

Two independent nodal reconstructions use the same original physical data. One uses the supplementary tensor-divergence evaluator with seven-point radial differences and Fourier angular differentiation. The second evaluates the explicit continuum Hessian formula above with analytic frame derivatives, direct seven-point radial second derivatives and Fourier angular derivatives. The routes share one-dimensional differentiation routines and are each independent of the Galerkin assembly. No manufactured solution is used. The measured physical residual is $\psi^{-10}(\operatorname{div}LW-8\pi G\psi^6S)$. Norms use the frozen $rq\psi^6$ weights with identical node sets for the two fields. Endpoints, first and last interior rows, two radial bands, bulk and full interior are reported. At 15 actual samples per grid, the closed projections and a separate direct 4D construction agree within $7.11\times10^{-15}$ for $E$, $1.12\times10^{-16}$ for $S_i$, $1.07\times10^{-14}$ for $S_{ij}$, and $2.14\times10^{-14}$ for its physical trace.

## Results on the original source

Values are weighted RMS from the independent explicit Hessian. Each cell gives archive → candidate.

| Region | 8 | 16 | 32 | 64 |
| --- | ---: | ---: | ---: | ---: |
| First interior row | 3.50e−3 → 1.99e−4 | 4.44e−3 → 8.81e−5 | 3.83e−3 → 3.20e−5 | 2.48e−3 → 1.16e−5 |
| Inner endpoint | 2.51e−2 → 7.06e−3 | 6.78e−2 → 2.74e−3 | 7.25e−2 → 9.50e−4 | 5.10e−2 → 3.33e−4 |
| Inner band | 3.20e−3 → 1.50e−4 | 4.43e−3 → 6.66e−5 | 3.73e−3 → 2.61e−5 | 2.34e−3 → 8.61e−6 |
| Bulk | 1.92e−3 → 8.58e−5 | 2.39e−3 → 1.39e−5 | 1.87e−3 → 3.33e−6 | 1.13e−3 → 8.43e−7 |
| Outer band | 9.95e−4 → 1.11e−4 | 7.67e−4 → 1.61e−5 | 5.55e−4 → 2.95e−6 | 3.28e−4 → 6.58e−7 |
| Outer endpoint | 1.72e−2 → 1.72e−4 | 1.09e−2 → 3.11e−5 | 4.04e−3 → 6.62e−6 | 1.08e−3 → 1.57e−6 |
| Whole interior | 1.94e−3 → 1.05e−4 | 2.39e−3 → 2.81e−5 | 1.92e−3 → 9.96e−6 | 1.17e−3 → 3.18e−6 |

Candidate Hessian residual orders are $1.17,1.46,1.46$ on the first row, $1.37,1.53,1.51$ at the inner endpoint, $1.17,1.35,1.60$ in the inner band, and $2.62,2.06,1.98$ in the bulk. The supplementary independent divergence route also improves on the archive, giving $5.79\times10^{-6}$ inner-band RMS at 64 compared with $1.07\times10^{-3}$. The roughly 210-fold reduction on the finest first row in the Hessian route is a strong sampled improvement. The remaining endpoint residual and limited row orders leave pointwise boundary accuracy unresolved.

To assess angular and mixed terms independently, compare the explicit continuum-action blocks at coincident nodes of successive grids. In the inner band the candidate’s mixed-block RMS differences are $7.42\times10^{-5},1.28\times10^{-5},2.03\times10^{-6}$, versus $7.76\times10^{-4},6.17\times10^{-4},2.33\times10^{-4}$ for the archive. Pure angular differences are $4.55\times10^{-5},3.82\times10^{-5},1.45\times10^{-5}$ versus $3.45\times10^{-5},3.07\times10^{-5},1.36\times10^{-5}$. There is no separate exact block target, so this tests refinement consistency rather than absolute error in each block. The candidate’s nested inner-band displacement differences fall from $5.95\times10^{-8}$ to $5.93\times10^{-9}$ and $4.84\times10^{-10}$. The archive gives $3.64\times10^{-7},2.03\times10^{-7},7.58\times10^{-8}$. These compare fields driven by each resolution’s own conformal factor.

| n | GMRES iterations | Actual unpreconditioned L2 | Relative L2 | Estimated Euclidean stiffness condition | Gauss-4/6 max displacement difference |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 8 | 1.30e−14 | 9.19e−11 | 87.7 | 7.26e−10 |
| 16 | 7 | 4.97e−14 | 6.81e−10 | 411 | 4.59e−12 |
| 32 | 7 | 1.09e−14 | 2.96e−10 | 1,824 | 2.83e−14 |
| 64 | 6 | 3.59e−14 | 1.94e−9 | 7,781 | 1.32e−16 |

All actual algebraic norms meet the implemented combined absolute/relative criterion. The absolute tolerance dominates in this weak load normalization. The eigenvalue estimates use the symmetric unpreconditioned stiffness in the Q2 nodal coefficient basis. The preconditioned spectral bound is analytical. The finest observed quadrature difference is within the scale where solve and floating-point errors can matter.

The unchanged gate was rerun as a separate process. The baseline is **12/12**. The diagnostic pairs with replaced momentum are **9/12**, with inner frozen momentum orders **0.4721, 0.8781, 1.1133**. The archived scalar samples and identical gate parameters were used. The gate’s Hamiltonian orders are essentially unchanged in this experiment, but its H residual includes the strain of the replaced W and is only a coupled-defect diagnostic here.

## Prior alternatives and what remains open

- The [first-row adjudication](../RUN43_FIRST_ROW_ADJUDICATION_20261103.md) and [product-rule obstruction](../RUN43_PRODUCT_RULE_OBSTRUCTION_20261105.md) establish the separate derivative-composition and angular multiplication defects. In the nodal algebra an idempotent $e_j$ would require $De_j=2e_jDe_j$, forcing every column of an exact Leibniz differentiation matrix to vanish. This pilot uses conforming test integration and applies no row matching or invented forcing.
- The [consistent momentum study](../RUN43_CONSISTENT_MOMENTUM.md) reached 12/12 with stabilized tensor divergence and scalar flux. Its fine physical inner-band residual by the supplementary metric is $1.07\times10^{-3}$. The [closure trials](../RUN43_COMPATIBLE_CLOSURE_20261103.md) reported 10/12 for compatible closure and 9/12 for second Hessian, coordinate dissipation, and the radial filter. Their archived fine supplementary inner residuals were $1.17\times10^{-3},8.56\times10^{-6},7.68\times10^{-5},7.63\times10^{-4}$, respectively. The [coordinated replay](../RUN43_COORDINATED_REPLAY_20261103.md) fourth-vector field already achieved $5.30\times10^{-7}$ with that metric and 9/12. Its finer residual is lower than this pilot’s. These historical solves used their own conformal factors and are contextual comparisons.
- The [exact mixed specification](../RUN43_CONSOLIDATED_MATHEMATICAL_SPECIFICATION_20261106.md) and [certified one-mesh assembly](../run43_assembly/README.md) established variational structure and linear stability. The [interval/Taylor pilot](https://github.com/GeraldDahlbergJr/dee-toroidal-helical-omniverse/blob/bbf804c9982bc20bb5dbd623d7d41142b0d5932c/research/run43_interval_taylor_pilot/REPORT.md) certified a Bernstein polynomial primitive and left the target error enclosure open. No manufactured-field evidence from those investigations enters the present physical validation.

The experiment isolates the momentum change by holding $\psi$ to its physical archive. In the full coupled problem $LW$ enters the Hamiltonian via $\psi^{-7}|LW|^2/8$, so a smaller momentum residual alone cannot establish coupled solution accuracy, positivity, or nonlinear convergence. The native Q2 strain is piecewise smooth with possible interelement jumps. Reconstructed nodal slopes cannot certify its pointwise strong residual or boundary traction. The original physics and data are preserved in the solved subsystem. A coupled follow-up is scientifically justified if an independent off-grid or weak residual assessment is combined with a scalar discretization and a positive coupled branch check under the original convergence requirements. The available evidence supports continued investigation of this operator and gives no basis for promotion.

## Reproduction and provenance

From repository root, on the research branch, with Python, NumPy and SciPy matching [environment.json](evidence/environment.json), use a fresh output directory.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
out=$(mktemp -d)
python research/run43_independent_momentum/experiment.py "$out/candidate"
python research/audit_run43_arrays.py "$out/candidate"  # expected exit 1, 9 of 12
archive="$PWD/research/results/coordinated-replay-20261103-1831/baseline"
mkdir "$out/baseline"
for n in 8 16 32 64; do
  ln -s "$archive/fields_$n.npz" "$archive/result_$n.json" "$out/baseline/"
done
python research/audit_run43_arrays.py "$out/baseline"    # expected exit 0, 12 of 12
```

The [evidence directory](evidence/) contains the new diagnostic arrays, hashes, environment, execution logs, both frozen audit outputs, comprehensive residuals and orders, and the nested comparisons. All archived input SHA-256 values were checked against their manifests. The experiment is isolated on the named branch from checkpoint `54183649dddf26f566b956d41618cec25b218212`. Only the new research directory is added.
