"""Read-only checkpoint/evidence verification; writes only the requested report."""
from pathlib import Path
import argparse, hashlib, json, platform, subprocess, sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
CHECKPOINT = '5ed79712e8170e1d794b358b16918094f79431f8'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path(__file__).with_name('checkpoint_verification.json'))
    args = ap.parse_args()
    tree = subprocess.check_output(['git', 'ls-tree', '-rz', CHECKPOINT], cwd=ROOT)
    files = {}; changed = []
    for entry in tree.split(b'\0'):
        if not entry: continue
        info, name = entry.split(b'\t', 1)
        mode, kind, blob = info.decode().split()
        rel = name.decode(); data = (ROOT/rel).read_bytes()
        actual = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        files[rel] = {'git_blob': blob, 'sha256': hashlib.sha256(data).hexdigest(), 'unchanged': actual == blob}
        if actual != blob: changed.append(rel)
    base = ROOT/'research/results/coordinated-replay-20261103-1831/baseline'
    archives = {}
    for n in (8,16,32,64):
        path = base/f'fields_{n}.npz'
        expected = json.loads((base/f'result_{n}.json').read_text())['sha256']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path) as a: u=a['psi']; w=a['W_cartesian']
        archives[str(n)] = dict(sha256=digest, expected_sha256=expected, matches=digest==expected,
            shape_valid=u.shape==(n+1,n,n) and w.shape==(3,n+1,n,n),
            finite=bool(np.isfinite(u).all() and np.isfinite(w).all()), psi_min=float(u.min()),
            psi_boundary_max_error=float(abs(u[[0,-1]]-1).max()), W_boundary_max_error=float(abs(w[:,[0,-1]]).max()))
    qdir=ROOT/'research/run43_nonpolynomial_validation'
    q=json.loads((qdir/'qualification.json').read_text())
    r=json.loads((qdir/'results.json').read_text())
    # Arithmetic consistency check of retained evidence, not a fixture rerun.
    fixture=[]
    for a,b,o in zip(q['rows'][1:],q['rows'][2:],q['orders']):
        ea=np.array(a['max_abs_error']); eb=np.array(b['max_abs_error'])
        ua=np.array(a['floating_error_bound']); ub=np.array(b['floating_error_bound'])
        low=np.log2((ea-ua)/(eb+ub))
        fixture.append(dict(pair=[a['n'],b['n']], lower_order=low.tolist(),
            recomputed_matches=bool(np.allclose(low,o['lower_order'],rtol=1e-14,atol=0)),
            qualified=bool(np.all(ea>10*ua) and np.all(eb>10*ub) and np.all(low>=1.5))))
    rates=[]
    for name in r['rows'][0]['regions']:
        for a,b in zip(r['rows'],r['rows'][1:]):
            x=a['regions'][name]; y=b['regions'][name]
            aa=np.array(x['rms']); bb=np.array(y['rms']); ua=np.array(x['uncertainty_rms']); ub=np.array(y['uncertainty_rms'])
            nominal=np.log2(aa/bb); valid=aa>ua
            low=[float(np.log2((aa[i]-ua[i])/(bb[i]+ub[i]))) if valid[i] else None for i in range(4)]
            stored=next(o for o in r['orders'] if o['region']==name and o['pair']==[a['n'],b['n']])
            matches=np.allclose(nominal,stored['nominal'],rtol=1e-14,atol=0) and all((v is None and w is None) or (v is not None and w is not None and np.isclose(v,w,rtol=1e-14,atol=0)) for v,w in zip(low,stored['conservative']))
            rates.append(dict(region=name,pair=[a['n'],b['n']],nominal=nominal.tolist(),sensitivity_adjusted=low,
                recomputed_matches=bool(matches), necessary_gate_passes=bool(all(v is not None and v>=1.5 for v in low) and np.all(ua/aa<=.1) and np.all(ub/bb<=.1))))
    manifest=json.loads((qdir/'input_sha256.json').read_text())
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in manifest.items()}
    report=dict(checkpoint=CHECKPOINT,python=sys.version,numpy=np.__version__,platform=platform.platform(),
        purpose='New checkpoint-wide preservation check and independent arithmetic verification of existing evidence; no solver or old evaluator execution.',
        tracked_file_count=len(files),all_preexisting_files_unchanged=not changed,changed=changed,
        input_hash_checks=hashes,accepted_archives=archives,fixture_evidence_arithmetic=fixture,
        baseline_evidence_arithmetic=rates,preexisting_files=files,new_DEE_solves=0,
        verdict='NOT YET ESTABLISHED',scope='Preservation and evidence arithmetic only; does not certify arbitrary-field interpolation error or physical accuracy.')
    passed=not changed and all(hashes.values()) and all(a['matches'] and a['shape_valid'] and a['finite'] and a['psi_min']>0 and a['psi_boundary_max_error']==a['W_boundary_max_error']==0 for a in archives.values()) and all(a['qualified'] and a['recomputed_matches'] for a in fixture) and all(a['recomputed_matches'] for a in rates)
    report['verification_status']='PASS' if passed else 'FAIL'
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('verification_status','checkpoint','tracked_file_count','all_preexisting_files_unchanged','new_DEE_solves','verdict')},indent=2))
    raise SystemExit(0 if passed else 1)
if __name__=='__main__': main()
