# Reproduce the Run 43 evidence checks

These commands evaluate archived evidence without running physical solves. From the root of a full clone, use Bash, Python 3.12.14, NumPy 2.3.5, jq, curl and GNU checksum tools. The replay used x86 extended precision with 63 explicit fraction bits in NumPy longdouble and one OpenBLAS thread. Keep the same shell for both blocks. The worktrees and output directories are disposable.

## Numerical replay

```bash
set -euo pipefail
git fetch origin
audit_work=$(mktemp -d)
git worktree add --detach "$audit_work/evidence" 19a0e7f94a16d759761a5fe55470263b837c9adf
git worktree add --detach "$audit_work/certification" 54183649dddf26f566b956d41618cec25b218212
export OPENBLAS_NUM_THREADS=1
cd "$audit_work/evidence"
replay=research/results/coordinated-replay-20261103-1831
sha256sum --quiet -c "$replay/SHA256SUMS"

mkdir "$audit_work/baseline" "$audit_work/coupled"
cp "$replay"/baseline/{fields_*.npz,result_*.json} "$audit_work/baseline/"
coupled=research/run43_coupled_recovery/results/fields
cp "$coupled"/{fields_*.npz,result_*.json} "$audit_work/coupled/"
python research/audit_run43_arrays.py "$audit_work/baseline"
cmp "$replay/baseline/independent-audit.json" "$audit_work/baseline/independent-audit.json"
cmp <(jq .checks "$replay/validation/baseline-independent-validation.json") \
    <(jq .checks "$audit_work/baseline/independent-audit.json")
coupled_status=0
python research/audit_run43_arrays.py "$audit_work/coupled" || coupled_status=$?
test "$coupled_status" -eq 1
cmp "$coupled/independent-audit.json" "$audit_work/coupled/independent-audit.json"

python research/run43_continuum_auditor_v2/validate.py "$audit_work/v2.json"
diff -u research/run43_continuum_auditor_v2/measurements.json "$audit_work/v2.json" || test "$?" -eq 1
```

The coupled exit code of 1 is expected and its output comparison must succeed. The v2 diff is inspected as data. On this replay the only differences were extended_residual_rms, float64_vs_extended_rms and evaluation_error_eps_scale_ratio on the coupled N=64 inner_band and inner_remaining rows. Relative variation reached about 1.21e−7. Classifications, scores, hashes and the other fields matched. Hybrid and 4N residual reconstruction columns rely on archived records.

## Hybrid archived RMS and manufactured certification accounting

Run the following from the evidence worktree left by the previous block. The certification scripts evaluate stored bounds and records.

```bash
python - <<'PY'
import json
import math
from pathlib import Path

root = Path("research")
hybrid = json.loads((root / "run43_hybrid/measurements.json").read_text())
rms = {row["n"]: row["band_rms"] for row in hybrid["experiments"]["hybrid"]["levels"]}
archived = json.loads((root / "run43_frozen_momentum/results/frozen-acceptance.json").read_text())
passed = 0
for check in archived["checks"]:
    coarse, fine = check["pair"]
    band, metric = check["band"], check["metric"]
    order = math.log2(rms[coarse][band][metric] / rms[fine][band][metric])
    accept = math.isfinite(order) and order >= 1.5
    assert abs(order - check["order"]) < 1e-12 and accept == check["passed"]
    passed += accept
    print(band, metric, coarse, fine, order, accept)
assert (passed, len(archived["checks"])) == (9, 12)
PY

cd "$audit_work/certification"
for dir in run43_bspline_completion run43_final_gate run43_final_certification; do
  (cd "research/$dir" && sha256sum --quiet -c SHA256SUMS)
done
study=research/run43_final_certification
python "$study/decide.py" > "$audit_work/decision.json"
cmp "$study/decision.json" "$audit_work/decision.json"
python "$study/nonlinear/certify.py" > "$audit_work/bounds.json"
cmp "$study/nonlinear/bounds.json" "$audit_work/bounds.json"
```

## Referee correction from the coupled evidence

```bash
cd "$audit_work/evidence"
python - <<'PY'
import csv
import math
from pathlib import Path

path = Path("research/run43_coupled_recovery/results/comparison.csv")
with path.open(newline="") as handle:
    rows = {int(row["n"]): row for row in csv.DictReader(handle)
            if row["field"] == "coupled" and row["region"] == "inner_band"}
primary = [float(rows[n]["continuum_rms"]) for n in (32, 64)]
alternate = [float(rows[n]["alternate_continuum_rms"]) for n in (32, 64)]
delta = float(rows[64]["reconstruction_delta_rms"])
assert abs(alternate[1] - primary[1]) <= delta <= alternate[1] + primary[1]
print(f"Finest scalar increase (%)     {100 * (alternate[1] / primary[1] - 1):.6f}")
print(f"Vector difference / primary (%) {100 * delta / primary[1]:.6f}")
print(f"Primary and alternate orders   {math.log2(primary[0] / primary[1]):.7f}, "
      f"{math.log2(alternate[0] / alternate[1]):.7f}")
PY
```

Expected results are 6.051653%, 60.870601% and orders 2.0661295 and 1.5248530. The reverse triangle inequality checks consistency of the distinct scalar and vector measurements.

## Zenodo deposit identity and checksums

Query the version records directly, use the provided content links and verify the deposited size and MD5 of each downloaded file. Review the date and SHA-256 outputs against the [inventory](REPORT.md).

```bash
cd "$audit_work"
python - <<'PY'
import hashlib
import json
from pathlib import Path
import subprocess

def download(url):
    return subprocess.check_output(["curl", "-fLsS", "--retry", "2", "--max-time", "60", url])


for record in ("24430854", "24330220"):
    metadata = json.loads(download(f"https://zenodo.org/api/records/{record}"))
    print(metadata["doi"], metadata["metadata"]["publication_date"])
    directory = Path(record)
    directory.mkdir(exist_ok=True)
    for entry in metadata["files"]:
        data = download(entry["links"]["self"])
        assert len(data) == entry["size"]
        assert "md5:" + hashlib.md5(data).hexdigest() == entry["checksum"]
        (directory / entry["key"]).write_bytes(data)
        print(entry["key"], hashlib.sha256(data).hexdigest())
PY
```
