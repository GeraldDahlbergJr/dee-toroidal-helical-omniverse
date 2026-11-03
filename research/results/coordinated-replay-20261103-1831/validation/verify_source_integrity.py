"""Compare original physical sources and audit with the preserved Nov 2 ref."""
from pathlib import Path
import subprocess,hashlib,json
r=Path('/workspace/dee-toroidal-helical-omniverse')
s=Path('/workspace/run43-investigation/solver/source-current')
out=Path(__file__).parent
paths=['docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md',
       'reproduce/dee_stress_energy_3p1.py',
       'reproduce/validate_toroidal_helical_dee_source.py',
       'reproduce/toroidal_helical_dee_source_run38_frozen_checkpoint.json',
       'reproduce/diagnose_toroidal_boundary_convergence.py']
rows=[]
for p in paths:
    raw=(s/p).read_bytes()
    old=subprocess.check_output(['git','-C',str(r),'show','origin/research/run43-consistent-momentum-20261102:'+p])
    blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    rows.append({'path':p,'sha256':hashlib.sha256(raw).hexdigest(),
                 'git_blob':blob,'identical_to_nov02_branch':raw==old})
assert all(x['identical_to_nov02_branch'] for x in rows)
assert rows[-1]['git_blob']=='5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
(out/'independent-original-source-integrity.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
