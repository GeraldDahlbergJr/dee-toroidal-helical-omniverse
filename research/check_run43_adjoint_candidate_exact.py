"""Bounded exact adjoint candidate certificate; no solver or physical-field imports."""
import json
from pathlib import Path
import sympy as s

R = s.Rational
n = 7
D = s.zeros(n)
D[0,0:3] = s.Matrix([[-R(3,2),2,-R(1,2)]])
D[-1,n-3:n] = s.Matrix([[R(1,2),-2,R(3,2)]])
for i in range(1,n-1):
    D[i,i-1], D[i,i+1] = -R(1,2), R(1,2)
H = s.diag(R(1,2),1,1,1,1,1,R(1,2))
I = s.eye(n)[:,1:-1]
B = s.diag(-1,0,0,0,0,0,1)
Div = H.inv()*(-D.T*H+B)
assert H*Div+D.T*H == B
K = I.T*H*Div*D*I
Gram = (D*I).T*H*(D*I)
assert K == -Gram and K == K.T
minors = [Gram[:j,:j].det() for j in range(1,n-1)]
assert all(x > 0 for x in minors)
# Exact local Taylor test on the actual smooth toroidal volume J=r(1+r*cos theta).
# Only rows 0,1,2,3 and nodes 0..4 enter these transpose rows; no grid sweep.
h = s.symbols('h', positive=True)
a,b,z = s.symbols('a b z', real=True)
J = lambda r: r*(1+r*z)
u = lambda r: (r-a)*(b-r)  # zero on BOTH original radial faces
up = lambda r: a+b-2*r
# Original first derivative is exact on this quadratic at all contributing rows.
F = [s.expand(J(a+j*h)*up(a+j*h)) for j in range(5)]
radial = [(-F[0]+F[2]/2)/(h*J(a+h)),
          (F[0]/4-F[1]/2+F[3]/2)/(h*J(a+2*h))]
leading = [s.simplify(s.limit(h*x,h,0)) for x in radial]
assert leading == [-(b-a)/2,(b-a)/4]
# Nonzero actual toroidal point theta=phi=0, c=(1,0,0): er=c,
# STF(er tensor c + c tensor er) contracted with er = (4/3)c.
actual = [s.simplify(x.subs({a:R(1,20),b:R(1,4),z:1})*R(4,3)) for x in leading]
assert actual == [-R(2,15),R(1,15)]
result = {
    'base_commit':'d206b2263bd0577c63026e1a3cc6ced5bb27a34d',
    'status':'STOP: strong boundary consistency disproven for stated candidate',
    'arithmetic':'exact SymPy rational/symbolic', 'sympy':s.__version__,
    'spatial_nodes':7, 'largest_assembled_matrix':[7,7],
    'normal_Green_identity':True, 'normal_Dirichlet_energy_identity':True,
    'normal_scalar_Gram_leading_minors':list(map(str,minors)),
    'full_toroidal_symmetry':'exact weighted algebraic identity in report',
    'full_toroidal_coercivity':'UNRESOLVED; normal reduction only certified',
    'boundary_h_times_radial_adjoint':list(map(str,leading)),
    'actual_inner_point_h_times_x_momentum_error':list(map(str,actual)),
    'boundary_accuracy':'first two interior rows have nonzero O(h_r^-1) defect',
    'weak_convergence':'not decided by this strong-consistency obstruction',
    'physical_solves':0, 'physical_field_comparison':'not reached: mandatory stop',
}
if __name__ == '__main__':
    print(json.dumps(result,indent=2))
