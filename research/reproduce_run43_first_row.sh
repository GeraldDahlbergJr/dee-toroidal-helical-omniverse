#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${DEE_ROW_TAG:?Set a fresh output tag}"
out="research/results/$DEE_ROW_TAG"
mkdir "$out"
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export DEE_SCALAR_ORDER=flux DEE_PRECONDITIONER=fourier
python research/certify_run43_first_row.py "$out/certificate" > "$out/certificate.log" 2>&1
python research/check_run43_radial_energy.py "$out/radial-energy.json" > "$out/radial-energy.log" 2>&1
python -u research/run43_manufactured_boundary_solve.py "$out/manufactured" > "$out/manufactured.log" 2>&1
python research/check_run43_physical_row_identity.py "$out/physical-row" research/results/radial-retry-20261103-fourth-vector > "$out/physical-row.log" 2>&1
python -m unittest discover -s tests -v > "$out/tests.log" 2>&1
find "$out" -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > "$out/SHA256SUMS"
