# Stage I: bounded exterior uniqueness investigation

2026-11-10 UTC. Exact Stage H parent
`72f2ff834bde7d046c5104310aa8b0b9df0254ee`;
branch `research/stage-i-exterior-spectral-uniqueness`.
One agent, sequential phases, no subagents or parallel tasks.

**Scientific classification: Outcome C. Exterior uniqueness is NOT
ESTABLISHED by this investigation.** The coefficient dictionary and two
local product identities are verified. The remaining proof obligation is
the exact parameter-uniform exterior estimate (I*) below, or a rigorously
certified imported version of it. This is a missing step in this audit,
not a claim that uniqueness is false, that no applicable theorem exists,
or that this is an unsolved problem in the mathematical literature.
No nonzero exterior solution or counterexample was constructed.

The companion `DEE_STAGE_I_SPECTRAL_THEOREM_AUDIT.md` records the actual
Koch–Tataru statements inspected, their asymptotic coefficient mapping,
source-sign inconsistency, retrieval limitations, and phase 1 decision.
This report does not silently correct the source's theorem sign or infer
an exterior statement from absence of global eigenvalues.

## The one selected method: cutoff exterior Carleman reduction

Set r=|x| and fix real Omega !=0. Use the exact H6 coefficients, the original
action and residual representation in the audit, without changing the
Einstein equations. Define

    Q=-div(A grad)-i Omega(2B.grad+div B)-Omega^2 C0.

Then Qz=0 is exactly the Stage H density equation multiplied by -1.
On a sufficiently far tail A is uniformly elliptic and real symmetric;
B,C0 are real, with the differentiated bounds given in the audit. For
z in H^1(E) intersect L^2(E), the equation is distributional and
z is H^2_loc there. Both norms use dx; the geometric measures are equivalent
on this end. No pointwise radiation condition, positive polynomial weight,
or superexponential decay is assumed.

Choose R>max(R0,1) sufficiently far into the elliptic H6 tail. Let chi be a
smooth real radial cutoff, zero for r<=R and one for r>=2R. Set v=chi z,
extended by zero through r<=R. This extension is H^1 and has no interface
jump because chi vanishes in a neighborhood of the interface. A smooth
elliptic continuation of A and real B,C0 inside the unused ball may be
chosen for notation only; it changes no physical field or physical equation.
No homogeneous equation is asserted for that extension of v inside its
source annulus. The exact commutator is

    Qv=F_R,
    F_R=-2 A^{ij}(partial_i chi)(partial_j z)
        -div(A grad chi)z-2i Omega(B.grad chi)z.                 (I1)

F_R belongs to L^2 and is supported in the fixed compact annulus
R<=r<=2R. Bounded coefficients and the original H^1 norm suffice. This
calculation includes the shift divergence and all angular tensor entries;
there is no radial reduction of z or A. Arbitrary finite H^{1/2} interface
trace of z at R0 is allowed. There is no boundary integral or zero-trace
assumption at that interface. Conversely, the existence of F_R is exactly
why a global eigenvalue exclusion cannot be applied directly to v.

The actual flat estimate (7) of Koch–Tataru Proposition 4 has an L^2 part
and an l^2 W^{1/4,4} part in three dimensions. For a constant positive
long-range V=Omega^2, A1 has tau0=0. Fix tau1=1/4. Their weights are

    h'_epsilon(t)=tau1+(tau exp(t/2)-tau1)tau^2/(tau^2+epsilon exp(t)),
    h_epsilon(0)=0,
    rho_epsilon(r)^4=h'_epsilon(log r)/r^2
       +h'_epsilon(log r)^2 (h''_epsilon(log r))_+/r^4.          (I2)

Here rho_epsilon is the **Carleman weight factor**, not the physical
amplitude rho in the Einstein–DEE action. For epsilon>0, the weight tapers
at infinity as in Proposition 4. Unweighted L^2 implies both
r^{tau1-1/2}v in L^2 and (1+r^2)^{tau1-1/2}v in L^2; the discrepancy between
the powers printed in Proposition 4 and Theorem 11 causes no integrability
obstacle for this choice. Proposition 4 itself explicitly accommodates
noncompactly supported v with its stated weighted condition. We do not
replace that condition by an unjustified limit of compact cutoffs with
exponentially growing errors.

The precise **unverified estimate needed for this route** is:

    || exp(h_epsilon(log r)) rho_epsilon v ||_{L^2}
       <= C || exp(h_epsilon(log r)) rho_epsilon^{-1} Qv ||_{L^2}.   (I*)

It must hold for the full real anisotropic H6 operator on a sufficiently far
tail, for v supported beyond a fixed radius with the above weighted
H^1_loc/L^2 admissibility and Qv a compactly supported L^2 source (in
particular v=chi z), uniformly for tau>=tau_* and
0<epsilon<=epsilon_*, with **one C independent of tau and epsilon**.
The fixed coefficient tail and support radius must be independent of those
parameters. If a compact-region weight modification is needed, its limiting
weight must retain a strict separation between the source annulus and every
farther fixed annulus, uniformly in tau. A precise admissibility/density
version is part of this obligation, not an implicit established estimate.
An L^p companion term may be retained, but is not needed for our compact
L^2 source once this L^2 inequality is available.

Section 7 explicitly points to an extension of (7), and the dictionary in
the audit checks (40)–(43). That is substantial evidence in favor of this
route. The displayed Theorem 11 sign, the global nature of Corollary 12,
and the unsupplied details of the modified exterior estimate prevent this
audit from certifying (I*) as stated. We neither assert that the extension
is false nor require a new theorem if an existing precise result suffices.

### What (I*) would prove, conditionally

This paragraph is an implication, **not a proof of (I*)**. For each fixed
tau, epsilon->0 gives

    h_0(log r)=2 tau(sqrt(r)-1),
    rho_0(r)^4=tau r^{-3/2}+(tau^3/2)r^{-5/2}.              (I3)

For v=chi z the right side of (I*) is on a fixed compact annulus; its
weight converges uniformly there. Fatou's lemma on the nonnegative left
side gives the corresponding limiting L^2 inequality. On any fixed
annulus K={a<r<b}, a>2R, the left weight is at least
c_K tau^{1/4} exp(2 tau(sqrt(a)-1)). On the source annulus,
rho_0^{-1}<=C_R tau^{-1/4} for tau>=1, so

    ||z||_{L^2(K)} <= C_{K,R} tau^{-1/2}
       exp[-2 tau(sqrt(a)-sqrt(2R))] ||F_R||_{L^2}.         (I4)

The constants may depend on the fixed background and annuli, not on tau.
Letting tau->infinity would imply z=0 on K. A countable annular exhaustion
would give zero exterior germ for r>2R. It would cover arbitrary interface
trace because (I1) used none. Extending that conclusion throughout a
connected uniformly elliptic E would require a valid local unique
continuation theorem with its hypotheses checked; it is not inferred
through an ergoregion or a degenerating principal form. No global interior
nonexistence claim is drawn even from this conditional implication.

I2–I4 are the explicit weight formula, its elementary limit, and an
inequality consequence under (I*), not new PDE estimates verified by a
symbolic program. The extra polynomial powers in I4 are conservative and
the exponential gap is positive since a>2R. This exhibits exactly what
uniformity is needed; a parameter-dependent remainder could destroy the
argument.

### Where the direct attempt stops

Attempting to deduce (I*) just by perturbing the flat Proposition 4 leaves

    div((A-I)grad v)

on its right side. Its second-order part is
(A-I)^{ij} partial_i partial_j v. Flat (7) controls v in a fractional
Sobolev norm and in the weighted L^2 norm, not this weighted second
derivative. The statement A-I=O(r^{-1}) therefore does not by itself
supply a norm-absorption inequality. Integration by parts against an
arbitrary weighted L^2 source norm does not remove that derivative loss.
This is a substantive failure of the direct flat-estimate perturbation,
not proof that the true variable-coefficient estimate fails. The latter
must treat A as the principal operator, which is precisely the additional
§7 step requiring certification.

For clarity, conjugation also leaves that principal tensor unchanged.
For a smooth real exponent phi, the exact identity is

    exp(phi) Q(exp(-phi)u)
      = Q u+2 A grad(phi).grad(u)
        +[div(A grad(phi))-grad(phi).A grad(phi)
                              +2i Omega B.grad(phi)]u.     (I5)

Thus a scalar weight does not turn the anisotropic second-order term into
a bounded multiplicative potential. The matched Hermitian shift is
retained. Smooth Coulomb terms may be accommodated in c in the coefficient
dictionary; no erroneous short-range L^2 classification is used to reject
them. This targeted attempt stops here. No alternative proof method,
scattering construction, numerical study, or symbolic search was pursued.

## Stopping decisions and scope of conclusions

1. Phase 1 read the actual relevant full theorem statements and proof
   sections, checked all coefficient and solution-space hypotheses, and
   found no fully certified direct exterior application in this audit.
2. Phase 2 made the one cutoff Carleman attempt above. It stopped at (I*),
   with an explicit derivative-loss reason why the flat estimate cannot
   simply be substituted. No repeated attempts followed.
3. Phase 3 consists only of the bounded symbolic verifier, unchanged
   historical regressions, provenance/preservation audit and publication.

**VERIFIED:** coefficient dictionary conditional on the stated H6 input;
cutoff and conjugation product identities; trace-free cutoff construction.
**CONDITIONAL:** (I*) would force a zero exterior germ as shown.
**NOT ESTABLISHED:** (I*), the full arbitrary-trace exterior uniqueness
conclusion, its extension to the compact interior, and H6 as a consequence
of the Einstein equations and finite energy alone.
**ERROR IDENTIFIED:** the opposite spectral sign in the inspected preprint's
(42); this is a source-applicability issue, not an Einstein–DEE counterexample
or a demonstrated error in Stage H. No error in the historical files was
repaired or any historical claim strengthened.

The Einstein equations give Stage H's leading constraints only conditional
on H6. Finite proper matter energy gives unweighted H^1 and L^2 for a fixed
nonzero frequency, not the uniform differentiated pointwise end expansions.
A separate weighted elliptic bootstrap remains necessary to justify H6 for
the full physical class. The exterior question and that bootstrap are
separate obligations. A proof only on a tail would also require additional
interior continuation hypotheses for global nonexistence. No complete
localized Einstein–DEE solution or independent peer review is claimed.

**Physical Review D readiness:** this is an auditable inconclusive technical
note. It does not support submission as a new rigorous general exterior
uniqueness theorem or a completed physical solution. Symbolic regression is
not peer review. Frozen recovery remains INCOMPLETE; physical tests NOT RUN.

## Verification and preservation record

The Stage I verifier checks exactly **two** local symbolic identities, I1
and I5, in three Cartesian dimensions with all six independent symmetric
principal entries and all three shift components as arbitrary functions.
It runs finite product differentiation and expansion only, writes no files,
and tests no global estimate, spectral theorem, numerical PDE or physical
model. H6 bounds, weighted inequalities, limit arguments and literature
applicability are analytical audit statements, not symbolic-test passes.

The records below are filled from the actual sequential command results
before commit. Publication itself does not require amending a report with
its own commit SHA; remote verification is returned in the final delivery.

### Actual sequential verification results

Python 3 with SymPy **1.14.0**. Stage I exited **0** and printed:

```text
SymPy 1.14.0
PASS: full anisotropic cutoff commutator
PASS: real exponential conjugation retains principal tensor
PASS: Stage I 2 exact local identities; exterior uniqueness NOT ESTABLISHED.
```

Its stderr was empty. The earlier verifiers and frozen package check were
then executed one at a time, with stdout, stderr and return codes captured
outside the repository. Every listed process completed; no test was skipped
or replaced by a historical pass report.

| Actual command | Exit code | Actual result |
| --- | --- | --- |
| `python DEE_STAGE_I_VERIFIER.py` | 0 | PASS, 2 local identities only |
| `python DEE_STAGE_H_VERIFIER.py` | 0 | PASS, 43 exact local checks |
| `python DEE_STAGE_G_VERIFIER.py` | 0 | PASS, 60 exact local checks |
| `python DEE_STAGE_F_VERIFIER.py` | 0 | PASS, 76 exact local checks |
| `python DEE_STAGE_E_VERIFIER.py` | 0 | PASS, 61 exact local checks |
| `python DEE_STAGE_D_VERIFIER.py` | 0 | PASS, 27 exact checks |
| `python DEE_STAGE_C_VERIFIER.py` | 0 | PASS, exact identities; analytical scopes unchanged |
| `python DEE_STAGE_B_VERIFIER.py` | 0 | PASS, exterior matching, currents, junction geometry, constraint and conservation |
| `python DEE_STAGE_A_HALF_VERIFIER.py` | 0 | PASS, all exact symbolic checks; no global existence/localization claim |
| `python DEE_STAGE_A_HALF_BOUNDED_VERIFIER.py` | 0 | PASS, bounded response, boundary and both phase-current identities |
| `python evidence/einstein_dee_sources/verify_evidence.py --package-only` | 0 | PASS, 49 package hashes, 20 archive members, pinned source copies |

All ten earlier-stage/package commands had empty stderr. Frozen stdout was
exactly:

```text
Git preservation NOT CHECKED in package-only mode
PASS: 49 package hashes, 20 archive members, pinned source copies
Recovery gate: INCOMPLETE. Physical tests: NOT RUN; approval required.
```

The historical approval note is quoted output; no physical experiment was
requested or run. **Frozen recovery remains INCOMPLETE. Physical tests
remain NOT RUN.** Provenance package-only mode does not verify Git-tree
preservation; the independent exhaustive check below does.

### Independent preservation and Git integrity

Before any additions, HEAD was verified as the exact Stage H parent:

    72f2ff834bde7d046c5104310aa8b0b9df0254ee
    tree 6107aab425810c45f9c263d12d90ba41ec2ae9f7

The baseline is **1,126 tracked paths**, counted directly from
`git ls-tree -r -z` rather than inferred from Stage H's report. Modes are
1,124 regular 100644 entries and two 100755 entries; no symlinks or gitlinks
are present in this baseline. The NUL-delimited recursive tree was saved
outside the repository before modifications and compared again against
the immutable parent. No AGENTS.md was found in the checkout or workspace
parent. The repository was cloned, then the authorized branch was created
directly at the specified parent; the different default-branch tip was
not used as a parent. The target remote branch was absent initially.

An independent external audit script enumerated **every** baseline path,
mode, type and object ID. It recomputed Git blob IDs from exact worktree
bytes (with symlink-target handling if such entries were encountered),
compared executable and symlink modes, and compared every historical entry
in the index and HEAD. It also walked every worktree file including ignored
files, excluding only .git, to detect unauthorized additions. Thus original
theory, all earlier stages/verifiers, all Run 43 evidence, frozen provenance
and every historical archive were covered without relying on filename
selection or package hashes.

Actual precommit audit results:

- `python /tmp/stage-i-preserve.py`: exit **0**. All **1,126** historical
  blobs and modes identical in worktree, index and HEAD; zero mismatches.
  **1,129** worktree files examined, including ignored files; exactly the
  three authorized new paths, all mode 100644.
- `git fsck --connectivity-only --no-dangling 72f2ff834bde7d046c5104310aa8b0b9df0254ee`:
  exit **0**, stdout and stderr empty. This checks object connectivity;
  exact file-content preservation is supplied by the independent blob audit.
- `git diff --check`: exit **0**, stdout and stderr empty.

The only additions are:

1. `DEE_STAGE_I_SPECTRAL_THEOREM_AUDIT.md`
2. `DEE_STAGE_I_EXTERIOR_UNIQUENESS.md`
3. `DEE_STAGE_I_VERIFIER.py`

No historical path is deleted, renamed, modified, chmodded, or replaced.
The index and completed commit are compared exhaustively again at publication;
any mismatch blocks commit or push. The final delivery reports the actual
commit SHA, remote branch/parent verification, three remote blob checks and
clean-worktree result, which cannot be embedded self-referentially in the
commit. No amend, rebase, force push, main merge, Stage J, manuscript work,
physical simulation or submission is authorized or performed. Stop after
Stage I publication.
