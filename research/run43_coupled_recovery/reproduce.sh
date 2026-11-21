#!/usr/bin/env bash
# Run from the repository root with the pinned packages in REPORT.md.
set -euo pipefail
work=${1:?provide a fresh work directory}
export OPENBLAS_NUM_THREADS=1
experiment=research/run43_coupled_recovery
mkdir -p "$work"
python "$experiment/manufactured.py" "$work/manufactured.json"
python "$experiment/solver.py" "$work/coupled"
status=0
python research/audit_run43_arrays.py "$work/coupled" || status=$?
# The gate uses 1 for a completed evaluation with failed checks.
if [[ $status -gt 1 || ! -s "$work/coupled/independent-audit.json" ]]; then
    exit 1
fi
research/run43_frozen_momentum/reproduce.sh "$work/prior" "$work/prior/replay"
python "$experiment/evaluate.py" "$work/coupled" "$work/prior" "$work/comparison.csv"
