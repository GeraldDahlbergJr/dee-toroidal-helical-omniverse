"""Verify every preexisting tracked blob and the saved physical input manifest."""
from pathlib import Path
import hashlib,json,subprocess,sys,platform
import numpy,sympy,mpmath
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
REF='5799e440b5f27bfe81c3ce5c44814a7d2d652157'
def git(*args):return subprocess.check_output(['git','-C',str(ROOT),*args])
entries=git('ls-tree','-r','-z',REF).split(b'\0')
checks={}
for entry in entries:
    if not entry:continue
    meta,path=entry.split(b'\t',1);mode,kind,blob=meta.split()
    path=path.decode();data=(ROOT/path).read_bytes()
    actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    checks[path]=actual==blob.decode()
manifest=json.loads((ROOT/'research/run43_final_validation/input_sha256.json').read_text())
inputs={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in manifest.items()}
archives={}
for n in (8,16,32,64):
    base=ROOT/'research/results/coordinated-replay-20261103-1831/baseline'
    expected=json.loads((base/f'result_{n}.json').read_text())['sha256']
    actual=hashlib.sha256((base/f'fields_{n}.npz').read_bytes()).hexdigest()
    archives[str(n)]={'sha256':actual,'matches_saved_result':actual==expected}
out={'reference_commit':REF,'branch':git('branch','--show-current').decode().strip(),'all_original_tracked_blobs_unchanged':all(checks.values()),'original_blob_count':len(checks),'changed_or_missing_originals':[p for p,v in checks.items() if not v],'original_input_manifest_count':len(inputs),'original_input_hash_checks':inputs,'accepted_archives':archives,'agent_count':1,'new_DEE_solves':0,'environment':{'python':platform.python_version(),'numpy':numpy.__version__,'sympy':sympy.__version__,'mpmath':mpmath.__version__},'evaluator_source_sha256':hashlib.sha256((ROOT/'research/run43_final_validation/evaluate.py').read_bytes()).hexdigest()}
assert out['all_original_tracked_blobs_unchanged'] and all(inputs.values()) and all(v['matches_saved_result'] for v in archives.values())
(OUT/'preservation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('original_input_hash_checks','accepted_archives')},indent=2))
