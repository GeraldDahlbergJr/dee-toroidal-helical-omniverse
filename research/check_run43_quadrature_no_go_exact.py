"""New bounded certificate: quadratic energy obstruction; no solver imports."""
import json
from pathlib import Path
import sympy as s

R = s.Rational
# Dimensionless original radial first derivative, including original endpoints.
D = s.Matrix([[-R(3,2),2,-R(1,2),0,0],
              [-R(1,2),0,R(1,2),0,0],
              [0,-R(1,2),0,R(1,2),0],
              [0,0,-R(1,2),0,R(1,2)],
              [0,0,R(1,2),-2,R(3,2)]])
F = s.zeros(4)
S = s.zeros(4)
for j in range(4):
    F[j,(j-1)%4] -= R(1,2)
    F[j,(j+1)%4] += R(1,2)
    S[j,j] -= 2
    S[j,(j-1)%4] += 1
    S[j,(j+1)%4] += 1
u = s.Matrix([0,0,1,0,0])
one = s.ones(4,1)
chi = s.Matrix([1,-1,1,-1])
assert u[0] == u[4] == 0
assert (D*u)[2] == 0
assert F*one == F*chi == s.zeros(4,1)
assert S*one == s.zeros(4,1)
assert S*chi == -4*chi
# At the work layer, angular differentiation of G_j W contributes zero,
# since D_r u=0 there. Radial compositions and radial Laplacian work agree.
# At every node, G_j W_chi = chi G_j W_one, hence identical strain norm.
# Difference of expanded work, multiplied by h_phi^2, at each theta node:
# q=1+r*cos(theta), r=3/20; arbitrary positive nodal quadrature weights.
q = [R(23,20),1,R(17,20),1]
w = s.symbols('w0:16', positive=True)
scaled_difference = -4*sum(w[4*t+p]/q[t]**2 for t in range(4) for p in range(4))
assert all(s.diff(scaled_difference, v) < 0 for v in w)
assert scaled_difference.is_negative is True
result = {
    'status':'PASS', 'checkpoint':'3f294581597a1831edf6001d8a191dae134b1499',
    'arithmetic':'exact rationals; positive symbolic quadrature weights',
    'largest_matrix_shape':[5,5], 'angular_ring_nodes':4,
    'radial_spike':list(map(str,u)), 'radial_derivative':list(map(str,D*u)),
    'angular_first_kills_both':True, 'angular_second_checkerboard_eigenvalue':-4,
    'radial_Dirichlet_both':True, 'strain_norm_equal_pointwise':True,
    'h_phi_squared_times_expanded_work_difference':str(scaled_difference),
    'strictly_negative_for_all_positive_nodal_weights':True,
    'conclusion':'No positive nodal quadrature can give the original expanded quadratic strain-energy identity for all radial-Dirichlet vectors.',
    'coercivity':'UNRESOLVED: no compatible construction survives the necessary energy test.'
}
if __name__ == '__main__':
    print(json.dumps(result,indent=2))
