"""Bitwise NPZ comparison after genuine fresh executions; no field rewriting."""
from pathlib import Path
import hashlib,json
s=Path('/workspace/run43-investigation/solver/source-current/research/results')
rows=[]
for variant,original in [('baseline','stabilized-flux-full'),
                         ('fourth','radial-retry-20261103-fourth-vector-repeat'),
                         ('local','boundary-local-20261103')]:
    for n in (8,16,32,64):
        a=s/f'live-20261103-{variant}'/f'fields_{n}.npz'
        b=s/original/f'fields_{n}.npz'
        da=hashlib.sha256(a.read_bytes()).hexdigest()
        db=hashlib.sha256(b.read_bytes()).hexdigest()
        rows.append(dict(variant=variant,n=n,fresh_sha256=da,
                         preserved_sha256=db,bitwise_identical=da==db,
                         preserved_path=str(b)))
output=Path(__file__).parent/'fresh-vs-preserved.json'
output.write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
