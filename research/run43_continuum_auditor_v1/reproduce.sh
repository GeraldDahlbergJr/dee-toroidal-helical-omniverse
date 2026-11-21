#!/usr/bin/env bash
# Run from the repository root with the versions in REPORT.md.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
work=$(realpath -m "${1:?provide a new scratch directory}")
study=research/run43_continuum_auditor_v1
base=3715241456898277e290fe135222582d538c182e
git diff --exit-code "$base" -- . ":(exclude)$study"
mkdir "$work"
research/run43_frozen_boundary/reproduce.sh "$work/frozen"
mkdir "$work/results"
python "$study/qualify.py" "$work/results/manufactured.csv"
python "$study/compare.py" "$work/frozen" "$work/results"
