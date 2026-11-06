# Reproduce the bounded diagnosis

Check out `run43-independent-physical-validation` from the repository and work at its root. The parent is `5ed79712e8170e1d794b358b16918094f79431f8`. No DEE solves are required or called. These commands read the accepted four baseline arrays and write only new evidence under `research/priority_one_validation/`.

Recorded software: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0. Use an isolated environment if dependencies are unavailable:

```sh
python -m venv /tmp/run43-validation-venv
. /tmp/run43-validation-venv/bin/activate
python -m pip install numpy==2.3.5 scipy==1.17.0 sympy==1.14.0
```

Verify preserved evidence first:

```sh
cd research/priority_one_validation
sha256sum -c SHA256SUMS
cd ../..
```

Run each bounded diagnostic once, sequentially:

```sh
OPENBLAS_NUM_THREADS=1 python research/priority_one_validation/agent1/diagnose_baseline.py > research/priority_one_validation/agent1/diagnostics.log 2>&1
OPENBLAS_NUM_THREADS=1 python research/priority_one_validation/agent2/decompose.py > research/priority_one_validation/agent2/decompose.log 2>&1
OPENBLAS_NUM_THREADS=1 python research/priority_one_validation/agent3/check_equations.py > research/priority_one_validation/agent3/equation_checks.log 2>&1
OPENBLAS_NUM_THREADS=1 python research/priority_one_validation/verify_checkpoint.py > research/priority_one_validation/checkpoint_verification.log 2>&1
```

Expected process exit code: 0 for all four. This indicates successful diagnostics, preservation and evidence arithmetic, **not** successful physical validation. Expected physical verdict: **NOT YET ESTABLISHED**. Preserved hashes and exact recorded results are in the committed JSON; floating-point diagnostics can vary slightly across platforms. The verification report includes runtime versions/platform, so its hash may differ on another system.

Agent 1 constructs derivative audit objects but calls no PDE solver. Its matched-closure subtraction evaluates the same immutable fields; it is not a corrected solve. Agent 2 reuses the preexisting standalone evaluator to isolate spread sources; this verifies decomposition, not a second wholly independent evaluator implementation. Agent 3 independently builds 4D stress projections and compares the production projection function. The main verifier reads all 952 starting tracked files, checks their Git blob identities, verifies NPZ/result hashes and exact nodal traces, and independently recomputes rate arithmetic from retained evidence. It does not rerun the historical fixture or claim new qualification.

Do not run historical renderers, solvers, `run43_final_validation/evaluate.py`, or the old qualification writers in place: they can overwrite preserved outputs. The historic nonpolynomial qualification remains inspectable at `research/run43_nonpolynomial_validation/`; its formulas, arithmetic enclosure and measured-grid limits are documented there. This run adds a diagnosis and an integrity/arithmetic verification, with a documented purpose for each new calculation.

To verify remote delivery, substitute the SHA from the final delivery message:

```sh
git rev-parse HEAD
git ls-remote origin refs/heads/run43-independent-physical-validation
```

Both must equal the delivered commit SHA. Stop after the diagnosis; no correction or new solve is implied by this reproduction procedure.
