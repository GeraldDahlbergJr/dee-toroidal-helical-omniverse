#!/usr/bin/env bash
# Run from the repository root with the dependencies in REPORT.md.
set -euo pipefail
work=${1:?provide a fresh work directory}
export OPENBLAS_NUM_THREADS=1
mkdir -p "$work"
study=research/run43_inner_asymptotic
prior=research/run43_frozen_momentum
coupled=research/run43_coupled_recovery/results/fields
"$prior/reproduce.sh" "$work/prior" "$work/prior/replay"
cp -a "$coupled" "$work/coupled-audit"
status=0
python research/audit_run43_arrays.py "$work/coupled-audit" || status=$?
[[ $status -eq 1 ]]
cmp "$coupled/independent-audit.json" "$work/coupled-audit/independent-audit.json"
python "$study/measure.py" "$work/prior" "$work/results"
