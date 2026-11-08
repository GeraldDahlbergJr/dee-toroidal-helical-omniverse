"""Extract original numbered TeX blocks; no algebra or physical tests."""
from pathlib import Path
import json,re,sys
root=Path(__file__).resolve().parent
source=root/'manuscripts/DEE_LaTeX_v2/main.tex'
text=source.read_text()
rows=[]
for m in re.finditer(r'\\begin\{equation\}.*?\\end\{equation\}',text,re.S):
    tag=re.search(r'\\tag\{([^}]+)\}',m.group())
    rows.append({'tag':tag.group(1) if tag else None,'line':text[:m.start()].count('\n')+1,'verbatim_tex':m.group()})
out=json.dumps({'source':'manuscripts/DEE_LaTeX_v2/main.tex','equations':rows},indent=2,ensure_ascii=False)+'\n'
dest=root/'EQUATIONS.json'
if '--check' in sys.argv:
    if dest.read_text()!=out: raise SystemExit('FAIL: equation extraction differs')
    print(f'PASS: {len(rows)} original equation blocks, no evaluation')
else: dest.write_text(out)
