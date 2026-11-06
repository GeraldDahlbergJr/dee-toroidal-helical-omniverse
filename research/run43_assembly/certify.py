"""One periodic brick: exact coefficient assembly and validated mass enclosure.

No physical load or physical solve. Existing integration evidence is read, never
regenerated. See README.md for the trace quotient and certificate transfer.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
from itertools import product
import sympy as s
from flint import arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / 'research/run43_integration'
spec = importlib.util.spec_from_file_location('validated_integration', OLD/'certify.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
ctx.prec = 128
u = s.symbols('u:3')
y = s.symbols('y:3')
xyz = s.symbols('x y z')
alphas = list(product(range(3), repeat=3))
mon = lambda a, variables: s.prod(v**k for v,k in zip(variables,a))
enrich = list(xyz) + [xyz[i]*xyz[j] for i in range(3) for j in range(i,3)]
E = [s.diag(1,-1,0), s.diag(1,1,-2)]
for i,j in [(0,1),(0,2),(1,2)]:
    e = s.zeros(3); e[i,j] = e[j,i] = 1; E.append(e)

def polynomial_coefficients(p):
    p = s.Poly(s.expand(p), *u)
    return s.Matrix([p.coeff_monomial(mon(a,u)) for a in alphas])

angular = [s.Integer(1), u[1]*(1-u[1])]
azimuthal = [s.Integer(1), u[2]*(1-u[2])]
full_q = [u[0]**a*angular[b]*azimuthal[c]
          for a,b,c in product(range(3),range(2),range(2))]
N = s.zeros(36,21)
for j,p in enumerate(full_q): N[:27,j] = polynomial_coefficients(p)
N[27:,12:] = s.eye(9)
I = s.zeros(21,4)
for j in range(4): I[4+j,j] = 1; I[8+j,j] = -1
J = s.eye(21)[:,list(range(8))+list(range(12,21))]
assert (I.row_join(J)).det() != 0

# Exact full strain coefficient map in the EXISTING small-cell tensor basis.
# xi=(y0+5)/20, eta=y1/(200*pi), zeta=y2/(200*pi).
Lraw = s.zeros(98,108)
substitution = dict(zip(u,[(y[0]+5)/20,y[1]/(200*s.pi),y[2]/(200*s.pi)]))
for j,a in enumerate(alphas):
    p = s.Poly(s.expand(mon(a,u).subs(substitution)), *y)
    assert s.expand(p.as_expr()-sum(p.coeff_monomial(mon(b,y))*mon(b,y)
                                   for b in alphas)) == 0
    for k,b in enumerate(alphas[1:]):
        c = p.coeff_monomial(mon(b,y))
        for i in range(3): Lraw[3*k+i,3*j+i] = c
for j,p in enumerate(enrich):
    for i in range(3):
        grad = s.Matrix([s.diff(p,x) for x in xyz])
        ei = s.eye(3)[:,i]
        T = ei*grad.T + grad*ei.T - s.Rational(2,3)*grad[i]*s.eye(3)
        coeff = [(T[0,0]-T[1,1])/2, -T[2,2]/2,T[0,1],T[0,2],T[1,2]]
        assert s.simplify(sum((c*e for c,e in zip(coeff,E)),s.zeros(3))-T) == s.zeros(3)
        for b,c in enumerate(coeff):
            q = s.Poly(c,*xyz)
            for k,v in enumerate([s.Integer(1)]+list(xyz)):
                Lraw[78+4*b+k,3*(27+j)+i] = q.coeff_monomial(v)
L = Lraw*s.kronecker_product(N,s.eye(3))
Ivec = s.kronecker_product(I,s.eye(3))
Jvec = s.kronecker_product(J,s.eye(3))
L0 = L*Ivec
assert L0.rank() == 12

tests = []
# Entire polynomial seam traces, not nodal comparisons.
for p in full_q:
    for k in [1,2]: assert s.expand(p.subs(u[k],0)-p.subs(u[k],1)) == 0
for j in range(4):
    p = sum(full_q[k]*I[k,j] for k in range(12))
    for t in [0,1]: assert s.expand(p.subs(u[0],t)) == 0
tests += ['entire_Q2_periodic_traces', 'entire_homogeneous_radial_traces',
          'interior_and_boundary_lift_direct_sum', 'strain_map_rank_12']
# All Cartesian quadratics are explicit raw generators (plus constant Q2).
for j,p in enumerate(enrich):
    for i in range(3):
        T = s.zeros(3)
        for b in range(5):
            q = sum(Lraw[78+4*b+k,3*(27+j)+i]*v
                    for k,v in enumerate([s.Integer(1)]+list(xyz)))
            T += E[b]*q
        g = s.Matrix([s.diff(p,x) for x in xyz]); ei = s.eye(3)[:,i]
        assert s.simplify(T-(ei*g.T+g*ei.T-s.Rational(2,3)*g[i]*s.eye(3))) == s.zeros(3)
assert Lraw[:,0:3] == s.zeros(98,3)
tests += ['all_30_Cartesian_P2_vector_strains_exact']
# Algebraic STF adjoint factor, independent of geometry and integration.
a = s.symbols('a:9'); b = s.symbols('b:9')
T = s.Matrix(3,3,a); G = s.Matrix(3,3,b)
T = (T+T.T)/2; T -= s.trace(T)*s.eye(3)/3
strain = G+G.T-s.Rational(2,3)*s.trace(G)*s.eye(3)
assert s.expand(s.trace(T.T*G)-s.trace(T.T*strain)/2) == 0
tests += ['STF_pointwise_adjoint_factor_one_half']

# Geometry/cofactor identity supplies the actual oriented face measure in the
# divergence theorem, including both paired chart seams.
rr,tt,pp = s.symbols('r theta phi', real=True)
X = s.Matrix([(1+rr*s.cos(tt))*s.cos(pp),
              (1+rr*s.cos(tt))*s.sin(pp),rr*s.sin(tt)])
for p in enrich:
    pullback = p.subs(dict(zip(xyz,X)), simultaneous=True)
    for angle in [tt,pp]:
        assert s.simplify(pullback.subs(angle,0)-pullback.subs(angle,2*s.pi)) == 0
# Nonconstant P2 restrictions on one radial torus are independent modulo
# constants. Laurent coefficients are exact full-face trace coefficients.
ztheta,zphi = s.symbols('ztheta zphi', nonzero=True)
cos_t,sin_t = (ztheta+1/ztheta)/2,(ztheta-1/ztheta)/(2*s.I)
cos_p,sin_p = (zphi+1/zphi)/2,(zphi-1/zphi)/(2*s.I)
faceX = [(1+s.Rational(1,20)*cos_t)*cos_p,
         (1+s.Rational(1,20)*cos_t)*sin_p,s.Rational(1,20)*sin_t]
traces = [s.Poly(s.expand(p.subs(dict(zip(xyz,faceX)),simultaneous=True)
                          *ztheta**2*zphi**2),ztheta,zphi) for p in enrich]
trace_coeff = s.Matrix([[p.coeff_monomial(ztheta**i*zphi**j) for p in traces]
                       for i,j in product(range(5),repeat=2) if (i,j)!=(2,2)])
assert trace_coeff.rank() == 9
tests += ['entire_Cartesian_enrichment_periodicity',
          'exact_radial_Fourier_trace_rank_9', 'all_Q2_affine_basis_changes_exact']
DX = X.jacobian([rr,tt,pp])
er = DX[:,0]; et = DX[:,1]/rr; ep = DX[:,2]/(1+rr*s.cos(tt))
for actual,target in [(DX[:,1].cross(DX[:,2]),-rr*(1+rr*s.cos(tt))*er),
                      (DX[:,2].cross(DX[:,0]),-(1+rr*s.cos(tt))*et),
                      (DX[:,0].cross(DX[:,1]),-rr*ep)]:
    assert (actual-target).applyfunc(s.trigsimp) == s.zeros(3,1)
tests += ['exact_oriented_cofactor_face_measures']

# Formal matrix identities hold for the same exact H, not for independently
# selected points of entry intervals. Hs is symmetric by construction.
Z = s.MatrixSymbol('Z',98,98); Hs = (Z+Z.T)/2
Ls = s.MatrixSymbol('L',98,63); Is = s.MatrixSymbol('I',63,12)
Bs = Ls*Is
assert (Bs.T*Hs*Bs/2).T == Bs.T*Hs*Bs/2
assert ((-Hs*Bs/2).T+Bs.T*Hs/2).expand() == s.ZeroMatrix(12,98)
assert (Is.T*Ls.T*Hs-Bs.T*Hs).expand() == s.ZeroMatrix(12,98)
assert L*Ivec == L0
tests += ['exact_C_symmetry', 'exact_variational_adjoint',
          'exact_reaction_lift_independence_modulo_interior_equation']

data = (OLD/'result.json').read_bytes()
assert hashlib.sha256(data).hexdigest() == 'f9b1e47b35608f4cc2028e0909ebe4d4c36caf437e488233c3e52b9f53a7cd57'
evidence = json.loads(data)
assert evidence['dimension'] == 98
assert evidence['certificate']['mass_eigenvalue_lower_power2'] == -332
tests += ['preserved_integration_evidence_hash_and_certificate']

def ball(expr):
    if expr == s.pi: return arb.pi()
    if expr.is_Rational: return arb(int(expr.p))/int(expr.q)
    if expr.is_Add: return sum(map(ball,expr.args),arb(0))
    if expr.is_Mul:
        v = arb(1)
        for x in expr.args: v *= ball(x)
        return v
    if expr.is_Pow: return ball(expr.base)**int(expr.exp)
    raise TypeError(expr)

# Bound the coefficient injectivity quantitatively with validated arithmetic.
B = arb_mat([[ball(L0[i,j]) for j in range(12)] for i in range(98)])
inverse = (B.transpose()*B).inv()
bound = sum((abs(inverse[i,j]) for i in range(12) for j in range(12)),arb(0))
k = 0
while not bound < arb(2)**k:
    k += 1
    assert k < 4096
tests += ['validated_L0_Gram_inverse_bound']

# Same validated integration procedure: one reference box, same entire tensor
# integrands and Frobenius weights. This is a coarse enclosure, never a rule
# substituted for the exact Lebesgue form. No historical integral is rerun.
box = [arb(1)/2+arb(0,arb(1)/2) for _ in range(3)]
small_coordinates = [20*box[0]-5,200*arb.pi()*box[1],200*arb.pi()*box[2]]
ts,_ = old.tensors(small_coordinates)
r = arb('1/20')+box[0]/5
q = 1+r*(2*arb.pi()*box[1]).cos()
assert r > 0 and q > 0
w = (4*arb.pi()**2/5)*r*q
H = [[arb(0) for _ in range(98)] for _ in range(98)]
for i in range(98):
    for j in range(i,98):
        v = w*sum((ts[i][c]*ts[j][c]*(2 if c>=3 else 1) for c in range(6)),arb(0))
        assert v.is_finite()
        H[i][j] = H[j][i] = v
HT = arb_mat(H)
LB = arb_mat([[ball(L[i,j]) for j in range(63)] for i in range(98)])
mixed = HT*LB
C = B.transpose()*HT*B/2
# Correlated exact definitions, rather than independent interval subtraction,
# preserve symmetry, adjointness and reactions even with coarse intervals.
for i in range(12):
    for j in range(i):
        assert C[i,j].overlaps(C[j,i])
        C[i,j] = C[j,i]
tests += ['finite_mass_enclosures', 'symmetric_mass_enclosures',
          'symmetric_condensed_enclosures']

def dyad(v): return list(map(int,v.man_exp()))
def encoded(M):
    return [[[dyad(M[i,j].mid()),dyad(M[i,j].rad())]
             for j in range(M.ncols())] for i in range(M.nrows())]
def exact(M): return [[str(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]

result = {
 'checkpoint':'6b71b0f3913454f1cff92f01233b1d5102b1bccc',
 'mesh':{'cells':1,'map':['1/20+xi/5','2*pi*eta','2*pi*zeta'],
         'periodic_face_pairs':[[1,0,1],[2,0,1]],'radial_zero_faces':[[0,0],[0,1]],
         'aspect_ratio':'10*pi < 32'},
 'dimensions':{'local_scalar':36,'periodic_scalar':21,'zero_trace_scalar':4,
               'tensor':98,'displacement':63,'zero_trace_displacement':12},
 'tensor_basis':'exact extension of run43_integration selected generators, same ordering',
 'raw_scalar_basis':[str(mon(a,u)) for a in alphas]+[str(p) for p in enrich],
 'periodic_scalar_basis':[str(p) for p in full_q]+[str(p) for p in enrich],
 'N':exact(N),'I':exact(I),'J':exact(J),'L':exact(L),'L0':exact(L0),
 'exact_definitions':{'H':'integral T_i:T_j dV','mixed':'H*L',
   'divergence_functional':'-L0.T*H/2','C':'L0.T*H*L0/2',
   'symmetric_mixed_block':'[[H,-H*L0],[-L0.T*H,0]]',
   'reaction':'Jvec.T*(f+L.T*H*sigma/2)',
   'reaction_lift_change':'Ivec.T*(f+L.T*H*sigma/2)=0 by interior equation'},
 'integration':{'precision_bits':128,'boxes':1,'method':evidence['method'],
    'entry_format':'[midpoint_dyadic,radius_dyadic], dyadic=[integer,exponent]',
    'uncertainty':'entrywise exact bounds; correlated algebraic definitions authoritative',
    'midpoint_positive_definiteness':'NOT CERTIFIED'},
 'H':encoded(HT),'mixed':encoded(mixed),'C':encoded(C),
 'certificate':{'preserved_result_sha256':hashlib.sha256(data).hexdigest(),
   'H_lower_power2':-332,'L0_Gram_inverse_norm_upper_power2':k,
   'C_lower_power2':-333-k,
   'method':'subdomain mass domination and validated coefficient injectivity',
   'scope':'exact H and exact C on this one mesh; symmetric mixed block is indefinite'},
 'tests':tests}
output = Path(__file__).with_name('result.json')
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'dimensions':result['dimensions'],'certificate':result['certificate'],
                  'tests':tests},indent=2))
