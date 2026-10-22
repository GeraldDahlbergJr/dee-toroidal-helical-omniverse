"""Run43 independent boundary validation; never promotes a failing candidate."""
import json, math, os, pathlib, subprocess, zipfile, hashlib, io
import numpy as np
from diagnose_toroidal_boundary_convergence import Audit
root=pathlib.Path(__file__).parent
artifact=int(os.environ["RUN43_ARTIFACT_ID"])
zip_path=root/"run43-final-stage.zip"
subprocess.run(["gh","api",f"repos/{os.environ['GITHUB_REPOSITORY']}/actions/artifacts/{artifact}/zip"],stdout=zip_path.open("wb"),check=True)
with zipfile.ZipFile(zip_path) as z:
    progress=json.loads(z.read("toroidal_direct_hessian_candidate_progress.json"))
    assert progress["status"]=="COMPLETE",progress.get("status")
    fine_bytes=z.read("toroidal_direct_hessian_fields_65_64.npz")
    fine=np.load(io.BytesIO(fine_bytes))
    psi,w=fine["psi"],fine["W_cartesian"]
    assert psi.shape==(65,64,64) and w.shape==(3,65,64,64)
    assert np.isfinite(psi).all() and np.isfinite(w).all() and np.min(psi)>0
    fine_row={"resolution":[65,64,64],"independent_audit":Audit(psi.shape).summarize(psi,w),"field_sha256":hashlib.sha256(fine_bytes).hexdigest()}
coarse=json.loads((root/"toroidal_direct_hessian_candidate.json").read_text())
rows=coarse["rows"]+[fine_row]
assert [r["resolution"][1] for r in rows]==[8,16,32,64]
regions=["inner_boundary_band","outer_boundary_band"]
metrics=["independent_flux_H","independent_tensor_M_order2"]
orders={}
gates={}
for region in regions:
    orders[region]={}
    for metric in metrics+["independent_tensor_M_order4"]:
        vals=[r["independent_audit"][region][metric]["volume_weighted_rms"] for r in rows]
        orders[region][metric]=[float(math.log2(a/b)) for a,b in zip(vals,vals[1:])]
    gates[region]=all(math.isfinite(p) and p>=1.5 for metric in metrics for p in orders[region][metric])
out={"status":"BOUNDARY_GATE_PASS_VALIDATION_PENDING" if all(gates.values()) else "BOUNDARY_GATE_NOT_MET","source_artifact_id":artifact,"source_run_id":38463009735,"source_progress_status":progress["status"],"criterion":"Independent H and second-order M >=1.5 on every refinement pair, in both boundary bands","rows":rows,"orders":orders,"boundary_gate":gates,"limitations":["Passing these tests is necessary but not sufficient for independent scientific validation.","Endpoint compatibility and full continuum verification remain separate."]}
path=root/"run43_boundary_validation_20261021.json"
path.write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({"status":out["status"],"orders":orders,"gate":gates},indent=2),flush=True)
