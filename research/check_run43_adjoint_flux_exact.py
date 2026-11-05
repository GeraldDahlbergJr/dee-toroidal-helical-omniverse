"""Exact four-node normal-stencil/Green-identity certificate; no solver imports."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = '4f8bcbe8e2ad2e33ce4ac100eb70b6d2c0b4ac73'
R = s.Rational
n = 4
D = s.Matrix([[-R(3,2),2,-R(1,2),0],[-R(1,2),0,R(1,2),0],
              [0,-R(1,2),0,R(1,2)],[0,R(1,2),-2,R(3,2)]])
H = s.diag(R(1,2),1,1,R(1,2))
Z = s.zeros(n)
G = [D,Z,Z]  # Cartesian constant-coefficient normal reduction only.
Hv = s.diag(H,H,H)
Ht = s.diag(*([H]*9))  # All nine tensor slots; Frobenius counts off-diagonals twice.
L = s.zeros(9*n,3*n)
Div = s.zeros(3*n,9*n)
B = s.zeros(3*n,9*n)
for i in range(3):
    for j in range(3):
        for k in range(3):
            block = (G[i] if j == k else Z) + (G[j] if i == k else Z)
            if i == j:
                block = block - R(2,3)*G[k]
            L[(3*i+j)*n:(3*i+j+1)*n,k*n:(k+1)*n] = block
        Div[i*n:(i+1)*n,(3*i+j)*n:(3*i+j+1)*n] = G[j]
    B[i*n,(3*i)*n] = -1
    B[(i+1)*n-1,(3*i+1)*n-1] = 1
# STF tensor injection: xx, yy, xy, xz, yz; zz=-xx-yy.
P = s.zeros(9*n,5*n)
for a in range(n):
    for slot,component,coefficient in [(0,0,1),(4,1,1),(8,0,-1),(8,1,-1),
                                      (1,2,1),(3,2,1),(2,3,1),(6,3,1),(5,4,1),(7,4,1)]:
        P[slot*n+a,component*n+a] = coefficient
assert s.ones(1,n)*D.T == s.zeros(1,n)
assert (Hv*Div + L.T*Ht/2-B)*P != s.zeros(3*n,5*n)
# Unique weighted adjoint with the stated B, on STF tensors.
Div_green = Hv.inv()*(-L.T*Ht/2+B)
assert (Hv*Div_green + L.T*Ht/2-B)*P == s.zeros(3*n,5*n)
J = s.zeros(3*n,3*(n-2))
for i in range(3):
    for a in range(1,n-1):
        J[i*n+a,i*(n-2)+a-1] = 1
K = J.T*Hv*Div_green*L*J
assert K == -(L*J).T*Ht*(L*J)/2
assert K == K.T
# Constant Cartesian test vectors give net flux, since L c=0.
C = s.diag(s.ones(n,1),s.ones(n,1),s.ones(n,1))
assert L*C == s.zeros(9*n,3)
assert (C.T*Hv*Div_green-C.T*B)*P == s.zeros(3,5*n)
# Original normal stencil has no positive diagonal SBP norm with endpoint-only flux:
# (H D + D.T H)[0,2] = -h0/2 since D[2,0]=0.
h0,h1,h2,h3 = s.symbols('h0 h1 h2 h3', positive=True)
Hs = s.diag(h0,h1,h2,h3)
normal_defect = (Hs*D+D.T*Hs)[0,2]
assert normal_defect == -h0/2
# Positive definite normal-reduction Dirichlet strain Gram matrix (Sylvester).
strain_gram = (L*J).T*Ht*(L*J)
minors = [strain_gram[:k,:k].det() for k in range(1,7)]
assert all(v > 0 for v in minors)
# Verify recovery, prior evidence hashes, and frozen source before emitting success.
assert subprocess.check_output(['git','rev-parse',CHECKPOINT],cwd=ROOT,text=True).strip() == CHECKPOINT
subprocess.run(['sha256sum','-c','research/results/product-rule-obstruction-20261105/SHA256SUMS'],cwd=ROOT,check=True,capture_output=True)
audit = 'reproduce/diagnose_toroidal_boundary_convergence.py'
blob = subprocess.check_output(['git','hash-object',audit],cwd=ROOT,text=True).strip()
assert blob == '5e4e54e1e1d84b785c4e66eac9b4d7fffdc43e22'
result = dict(status='PASS',arithmetic='exact SymPy rationals and symbolic positive weights',
    sympy=s.__version__,nodes=n,checkpoint=CHECKPOINT,frozen_audit_blob=blob,
    scope='Constant-coefficient normal reduction plus general symbolic Green derivation; no toroidal-grid stability or physical validation',
    D=[[str(x) for x in D.row(i)] for i in range(n)],
    original_normal_SBP_obstruction=str(normal_defect),
    weighted_Green_STF_identity=True,Dirichlet_energy_identity=True,
    constant_test_flux_balance=True,normal_Dirichlet_Gram_leading_minors=[str(v) for v in minors],
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
out = Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=False)
(out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
