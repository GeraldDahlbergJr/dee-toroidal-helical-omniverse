# Run 43 continuity handoff

Evidence checkpoint: 7325338b2824547bb34b7da875caef541a1e77fd
Branch: run43-independent-physical-validation
Decision: PROCEED to a bounded qualification experiment after authorization.
Independent physical validation remains NOT YET ESTABLISHED. PRD submission
readiness is on hold. No experiment, implementation, or manuscript is authorized
by this handoff.

## Verified evidence

- Frozen baseline and fresh reproduction: 12/12. Nonpolynomial evaluator fixture:
  PASS. Independent continuum validation: NOT YET ESTABLISHED.
- Finest inner-band momentum rates: 1.027781, 1.027781, 0.990610. Radial
  sensitivity approximately 74–79%; finest inner H sensitivity is also 26.8%.
- First interior radial composition defect: h f'''/4, with opposite sign at the
  other face, multiplied in momentum by I+er er^T/3. The frozen audit shares
  this closure. At n=64 the two inner-band terms are each ~74 times the frozen
  residual, with cosine −.99990836. The frozen residual tracks angular
  stabilization to solve error.
- Source and stress-energy algebra checked. No established source correction.
- Preserve all historical evidence. Do not rerun matched endpoints, local/global
  radial replacements, direct Hessians, fourth-vector, compatible closures,
  coordinate KO, or radial filters. The one-brick mixed certificate does not
  qualify a physical solve.

## One candidate

Use a primal conforming Galerkin discretization of both constraints with
degree-four tensor B-splines in the exact toroidal chart, open radial knots,
simple interior knots and periodic angular bases. Use Cartesian vector
components and impose psi=1 and W=0 strongly on both radial faces.
Use the Euclidean chart volume rq and exact chart gradient. No stabilization,
extra boundary condition, or change to the physical source.

Test -Delta and -div L with the scalar gradient form and the vector form
a(U,V)=integral LU:LV/2. For radial-zero periodic fields, integration by
parts gives a(U,U)=||grad U||²+||div U||²/3. Thus the exact momentum block
is coercive uniformly in mesh size, as is the scalar principal block.
Galerkin orthogonality and conforming approximation give consistency for
these blocks. This is neither a coupled nonlinear stability proof nor a
prediction of passing strong-residual rates. Verify quadrature and local
coupled conditioning before any claim about solution error.

Use Delta psi + psi^(-7)|LW|²/8 + 2 pi G E psi^5=0 and
div LW=8 pi G psi^6 S. Retain G=.001, R=1, r=.05 to .25,
q=1+r cos(theta), E=.0127638+.068614375*(r^-2+9q^-2)*psi^-4,
S=-.030316*(e_theta/r-3e_phi/q). Pin exact source inputs to the archived
implementation, with the displayed decimals as cross-checks.

## Next action after authorization

One isolated implementation, one nonpolynomial manufactured coupled fixture,
n=8,16,32. No physical DEE solves in this first authorization. Require correct
source conventions, exact trace/seam construction, no unexpected kernel,
relative linearized algebraic residual <=1e-10, and H and all Cartesian M
rates >=1.5 on both pairs in the inner band, bulk, outer band and first
interior planes. Score trace defects <=100 machine eps times declared field
scales. Use exact analytic fields for fixture field-error checks. Stop on a
miss or budget exhaustion. Suggested cap: one CPU thread, 30 minutes total,
4 GiB, no dense assembly or retries.

Use independently coded Cartesian chain-rule differentiation of saved spline
coefficients and independent residual integration. Document exact spline
representation, compare with independent derivative evaluation and quadrature,
and justify bounds on evaluation/integration and solver effects. The 10%
budget applies to each scored physical metric; diagnostic agreement alone
is not an error bound. Use p_low=log2((Ec-uc)/(Ef+uf)), positive numerator
required. Freeze dimensional scales and accuracy tolerances before scoring.
A spline defines an approximate field without a choice of nodal reconstruction;
this does not bound error relative to an unknown continuum solution.

After fixture qualification, request separate authorization for original-source
grids 8,16,32, with 64 only if needed and separately budgeted. Stop at first
failed gate. Require positivity, radial and angular conditions, conservative
orders >=1.5, uncertainty <=10%, and predeclared absolute residual targets
including first interior planes. Report the frozen audit separately. Require
a verified local inverse/nonlinearity estimate or comparable a posteriori
result for any absolute solution-error claim. Independently reproduce any
successful physical result from the pinned commit before manuscript work.

## Reproduction and evidence

From the pinned checkout, verify sealed evidence:
    (cd research/priority_one_validation && sha256sum -c SHA256SUMS)

See priority_one_validation/REPRODUCE.md (beside this file) for Python versions
and the four no-solve diagnostic commands. Use a disposable checkout for any
replay because those scripts write diagnostic outputs. See REPORT.md and
agent1/REPORT.md through agent3/REPORT.md for the diagnosis, radial sensitivity
and equations. Historical evidence is indexed in agent1/REPORT.md and in
research/RUN43_COORDINATED_REPLAY_20261103.md and
research/RUN43_COMPATIBLE_CLOSURE_20261103.md. The qualified evaluator evidence
is in research/run43_nonpolynomial_validation/.
