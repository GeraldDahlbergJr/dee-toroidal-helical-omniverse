"""Tiny angular algebra certificate; no fields, physical forcing, or solves."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

root = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=False)
n = 16
h = 2*np.pi/n
x = np.arange(n)*h
eye = np.eye(n)
D = (np.roll(eye, -1, axis=0)-np.roll(eye, 1, axis=0))/(2*h)
A = (np.roll(eye, -1, axis=0)+np.roll(eye, 1, axis=0))/2
# Basis pairs certify the bilinear identity for every pair of nodal arrays.
ordinary = shifted = 0.
for a in eye:
    for f in eye:
        lhs = D@(a*f)
        ordinary = max(ordinary, float(np.max(abs(lhs-(D@a)*f-a*(D@f)))))
        shifted = max(shifted, float(np.max(abs(lhs-(D@a)*(A@f)-(A@a)*(D@f)))))
assert shifted < 1e-13 and ordinary > 1
k, m = 1, 2
a, f = np.exp(1j*k*x), np.exp(1j*m*x)
defect = D@(a*f)-(D@a)*f-a*(D@f)
symbol = 1j*(np.sin((k+m)*h)-np.sin(k*h)-np.sin(m*h))/h
symbol_error = float(np.max(abs(defect-symbol*a*f)))
assert symbol_error < 1e-13
audit = 'reproduce/diagnose_toroidal_boundary_convergence.py'
blob = subprocess.check_output(['git', 'hash-object', audit], cwd=root, text=True).strip()
assert blob == '5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
result = dict(scope='Finite angular algebra only; no physical convergence result',
    checkpoint='ff07102de88a4b63686a4f71f9fa851602ef7e0b', frozen_audit_blob=blob,
    angular_nodes=n, spacing=h, basis_pairs=n*n,
    ordinary_product_rule_max_defect=ordinary,
    shifted_product_rule_max_defect=shifted,
    mode_pair=[k,m], mode_defect_magnitude=float(abs(symbol)),
    mode_identity_max_error=symbol_error, numpy=np.__version__,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'certificate.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
