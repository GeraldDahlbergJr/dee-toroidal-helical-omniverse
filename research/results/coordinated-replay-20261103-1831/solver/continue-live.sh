#!/usr/bin/env bash
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier
# The baseline process was launched separately; wait for its captured exit.
while [[ ! -f live-evidence/baseline-solve.exit ]]; do sleep 5; done
for variant in fourth local; do
  if [[ "$variant" == fourth ]]; then
    script=research/run43_fourth_vector.py
  else
    script=research/run43_local_boundary.py
  fi
  output="research/results/live-20261103-$variant"
  if [[ -e "$output" ]]; then
    echo "Refusing to overwrite $output" >&2
    exit 1
  fi
  DEE_EXPERIMENT="live-20261103-$variant" python -u "$script" > "live-evidence/$variant-solve.log" 2>&1
  status=$?
  echo "$status" > "live-evidence/$variant-solve.exit"
done
