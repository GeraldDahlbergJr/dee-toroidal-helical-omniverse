# Run 43 continuation: angular product-rule obstruction

No new physical candidate or full-grid run is justified. This continuation starts
at coordinated checkpoint `ff07102de88a4b63686a4f71f9fa851602ef7e0b`, on separate
branch `research/run43-product-rule-obstruction-20261105`. No old experiment was
rerun. Read: coordinated replay report, its operator and boundary reports,
finite-grid-terms report, endpoint budget, and compatible-closure report.

## New analytical result applied to this problem

An exact ordinary Leibniz rule cannot be restored by choosing a better linear
nodal differentiation matrix while retaining all pointwise products. This is a
finite-dimensional algebra obstruction, not a new theorem about continuum
derivatives. If D(af)=(Da)f+a(Df) for all nodal arrays, take a=f=e_j, a nodal
idempotent. Since e_j squared equals e_j, De_j=2 e_j De_j. Off node j this
forces De_j=0; at j it also forces De_j=0. Every column vanishes, so D=0.
The result applies to finite Fourier matrices as well as finite differences;
spectral differentiation only obeys the ordinary rule on suitable unaliased
subspaces. It does not exclude tensor identities holding on particular fields
or cancellation between terms.

There is instead an exact centered angular identity. Write

    D=(Splus-Sminus)/(2h), A=(Splus+Sminus)/2.
    D(af)=(Da)(Af)+(Aa)(Df).
    [D,M_a]-M_Da = M_(Aa-a) D + M_Da (A-I).

M means pointwise multiplication. This follows by expanding the two shifts,
without Taylor assumptions. Thus the discrepancy includes averages of both the
coefficient and differentiated field. For smooth functions its leading ordinary
product-rule defect is h^2*(a'' f'+a' f'')/2. For Fourier modes k,m its exact
symbol is i*(sin((k+m)h)-sin(kh)-sin(mh))/h. It is nonzero even without aliasing.
The toroidal basis and reciprocal scale factors are variable multiplication
coefficients, so demanding an ordinary discrete chain/product rule for all
Cartesian component arrays cannot remove their commutators.

This sharpens the November 3 diagnosis rather than redoing its measured error
budget. That budget already established material radial, angular, and mixed
composition terms. Radial endpoint repair alone cannot restore an ordinary
product rule globally. A valid mimetic method would need compatible averages,
products, and fluxes, plus independent consistency and stability justification.
The shifted identity alone does not select such a method or justify inserting
the frozen operator defect into an accurate physical solve. It does not prove
that the required twelve-check gate is impossible.

## Continuum source compatibility

For the stated smooth periodic shell with W=0 on both radial boundaries,
integration by parts gives

    integral W dot div(LW) dV = -1/2 integral |LW|^2 dV.

The boundary term vanishes because W vanishes, not because traction is prescribed.
The Dirichlet longitudinal operator has no nonzero conformal Killing field in
its kernel on this shell; the usual Korn/coercivity argument gives a unique weak
momentum solution for prescribed admissible psi and source. Therefore a nonzero
net matter momentum is not itself a source incompatibility: boundary traction
can carry its integral. Imposing zero traction or W'''=0 as an additional
condition would change this problem. This standard continuum argument does not
prove solvability of the nonlinear Hamiltonian coupling, finite-grid stability,
or convergence of any candidate. No such claims are made here.

## Focused certificate and decision

`check_run43_product_rule.py` uses one 16-node angular ring, not radial fields or
synthetic physical solutions. All 256 nodal basis pairs certify the bilinear
shifted identity: maximum defect 0. The ordinary identity defect is
1.2732395447351628. The unaliased mode pair k=1,m=2 has defect magnitude
.4224880799890335; the exact symbol agrees within 1.50e-15. These are algebra
diagnostics, not H/M2 norms or convergence results. The script also verifies
that the frozen audit blob is unchanged.

No solver correction was derived with enough accuracy/stability evidence to
warrant even a new physical solve. No n=8,16,32,64 sweep was launched. Existing
baseline 12/12, fourth 9/12, and local 8/12 outcomes remain unchanged; the other
November 3 compatible-closure failures were read, not repeated. All twelve
orders >=1.5 in both original physical bands remain mandatory.

The next blocker is a whole-system boundary/flux/product discretization with an
independent error and stability justification and a focused diagnostic showing
improvement on the original fields. Exact ordinary nodal Leibniz restoration is
now excluded as a general construction. This result supplies no new physical
validation, and initial-data constraints do not establish Einstein-DEE evolution.

## Preservation

Only this report, the new certificate script, and its JSON/hash manifest are
added. Original sources, frozen audit, checkpoints, and historical results are
untouched. Evidence is small and preserved through Git on the separate research
branch, the same repository-based method used by the checkpoint. GitHub Actions
artifact storage was not tested or assumed available; no expensive job needs it.
The old local Git bundle is not present in this fresh workspace; it was not
claimed recovered. The remote checkpoint and its tracked evidence were cloned.

Reproduce with a fresh output path:

    python research/check_run43_product_rule.py research/results/UNIQUE-product-rule

The script refuses to overwrite an existing directory.
