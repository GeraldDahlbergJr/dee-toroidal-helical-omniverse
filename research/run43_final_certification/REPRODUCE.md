# Reproduce without PDE solves

Use Python 3.12, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0. Start at the
repository root. The commands below perform checks and archived-field
evaluations only. Fresh empirical runs contain different timing fields;
compare scores and hashes, not elapsed time. No installation was needed.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
(cd research/run43_bspline_completion && sha256sum -c SHA256SUMS)
(cd research/run43_final_gate && sha256sum -c SHA256SUMS)
(cd research/run43_final_certification && sha256sum -c SHA256SUMS)
python research/run43_final_certification/build_budget.py > /tmp/run43-budget.json
cmp /tmp/run43-budget.json research/run43_final_certification/error_budget.json
python research/run43_final_certification/decide.py > /tmp/run43-decision.json
cmp /tmp/run43-decision.json research/run43_final_certification/decision.json
python research/run43_final_certification/nonlinear/certify.py > /tmp/run43-bounds.json
cmp /tmp/run43-bounds.json research/run43_final_certification/nonlinear/bounds.json
python research/run43_final_certification/audit/verify.py > /tmp/run43-independent.json
cmp /tmp/run43-independent.json research/run43_final_certification/audit/independent.json
python research/run43_final_certification/audit/denominator.py > /tmp/run43-denominator.json
cmp /tmp/run43-denominator.json research/run43_final_certification/audit/denominator.json
python research/run43_final_certification/residual/pilot.py
python research/run43_final_certification/residual/coarse_pilot.py
CERT_RULE=6 timeout --kill-after=5s 240 python research/run43_final_certification/residual/reproduce.py
CERT_RULE=8 timeout --kill-after=5s 240 python research/run43_final_certification/residual/reproduce.py
timeout --kill-after=5s 240 python research/run43_final_certification/audit/replay.py > /tmp/run43-replay.json
cmp /tmp/run43-replay.json research/run43_final_certification/audit/replay.json
timeout --kill-after=5s 240 python research/run43_final_certification/nonlinear/weak_defect.py 32
timeout --kill-after=5s 240 python research/run43_final_certification/nonlinear/weak_defect.py 32 8
timeout --kill-after=5s 240 python research/run43_final_certification/nonlinear/weak_defect.py 32 6 --compare
```

The flag `--compare` reassembles both rules without solving. It is an
empirical assembly diagnostic. Replays overwrite their workstream JSON with
fresh timing metadata; verify the evidence seal before replay, or use a
disposable copy. Historical residual q6/q8 cumulative elapsed times were
about 17.4/36.0 seconds across all grids; diagnostic timings depend on hardware.

Confirm baseline preservation with:

```sh
git diff --exit-code c5580eea5b216052ba5e91c6725e0304a7852a99 -- research/run43_bspline_completion research/run43_final_gate research/run43_bspline_qualification
```
