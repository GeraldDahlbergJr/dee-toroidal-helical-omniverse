"""Read-only provenance verifier. Never import or execute scientific code."""
from pathlib import Path
import hashlib,json,subprocess,sys,zipfile
P=Path(__file__).resolve().parent
def sha(b): return hashlib.sha256(b).hexdigest()
def blob(b): return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def require(ok,message):
    if not ok: raise SystemExit('FAIL: '+message)
count=0
for line in (P/'SHA256SUMS').read_text().splitlines():
    digest,name=line.split('  ',1)
    require(sha((P/name).read_bytes())==digest,name);count+=1
inventory=json.loads((P/'PINNED_TREE.json').read_text())
files=[e for e in inventory['entries'] if e['type']=='blob']
for e in files:
    q=P/'repository'/e['path']
    if q.exists(): require(blob(q.read_bytes())==e['sha'],'source copy '+e['path'])
with zipfile.ZipFile(P/'manuscripts/DEE_LaTeX_v2(5).zip') as z:
    members=[i for i in z.infolist() if not i.is_dir()]
    for i in members: require((P/'manuscripts'/i.filename).read_bytes()==z.read(i),'archive '+i.filename)
if '--package-only' not in sys.argv:
    root=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=P,text=True).strip())
    require(subprocess.check_output(['git','rev-parse',inventory['commit']],cwd=root,text=True).strip()==inventory['commit'],'pinned commit')
    subprocess.run(['git','merge-base','--is-ancestor',inventory['commit'],'HEAD'],cwd=root,check=True)
    raw=subprocess.check_output(['git','ls-tree','-rz','HEAD'],cwd=root)
    head={}
    for row in raw.split(b'\0'):
        if row:
            meta,path=row.split(b'\t',1);mode,typ,oid=meta.decode().split();head[path.decode()]=(mode,oid)
    for e in files:
        require(head.get(e['path'])==(e['mode'],e['sha']),'HEAD changed '+e['path'])
        require(blob((root/e['path']).read_bytes())==e['sha'],'working file changed '+e['path'])
    prefix=P.relative_to(root).as_posix()+'/'
    require(all(n in {e['path'] for e in files} or n.startswith(prefix) for n in head),'addition outside evidence package')
    print(f'PASS: all {len(files)} pinned files unchanged in HEAD and working tree')
else: print('Git preservation NOT CHECKED in package-only mode')
print(f'PASS: {count} package hashes, {len(members)} archive members, pinned source copies')
print('Recovery gate: INCOMPLETE. Physical tests: NOT RUN; approval required.')
