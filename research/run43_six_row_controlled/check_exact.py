"""One six-row family; exact small matrices only. No solver imports."""
import sympy as s
import json
from pathlib import Path
n=16; m=6
H=s.diag(s.Rational(1,2), *([1]*(n-2)), s.Rational(1,2))
x=s.Matrix(range(n)); one=s.ones(n,1)
D=s.zeros(n); D[0,0]=-1
q={}
for i in range(m):
 for j in range(i+1,m):
  q[i,j]=s.Symbol(f'q{i}{j}')
  D[i,j]=q[i,j]/H[i,i]; D[j,i]=-q[i,j]/H[j,j]
for i in range(m,n-m):
 D[i,i-1]=-s.Rational(1,2); D[i,i+1]=s.Rational(1,2)
D[m-1,m]=s.Rational(1,2)
# Outer block is reflected, including its centered transition coupling.
for i in range(m):
 for j in range(m+1): D[n-1-i,n-1-j]=-D[i,j]
eq=list(D*one)+list(D*x-one)
sol=s.linsolve(eq,list(q.values())); point=next(iter(sol))
D=D.subs(dict(zip(q.values(),point)))
assert D*one==s.zeros(n,1) and D*x==one
B=s.diag(-1,*([0]*(n-2)),1)
assert H*D+D.T*H==B
x2=x.applyfunc(lambda z:z*z)
e=s.simplify(D*x2-2*x)
eL=s.Matrix([e[i] if i<m else 0 for i in range(n)])
r=s.simplify(D*eL)
assert all(r[i]==0 for i in range(m+1,n))
# All inner affected/transition rows 1..6 must vanish.
conditions=[s.factor(r[i]) for i in range(1,m+1)]
invariant=s.factor(sum(i*H[i,i]*r[i] for i in range(1,m+1)))
assert invariant==-s.Rational(1,2)
assert s.factor((one.T*H*eL)[0])==s.Rational(1,2)
# A single member witnesses that the old kernel is not universal.
free=sorted(D.free_symbols,key=str)
D0=D.subs(dict.fromkeys(free,0))
I=s.eye(n)[:,1:n-1]
rank=(D0*I).rank()
old=s.zeros(n,1); old[1]=1; old[2]=-1
out={'base':'80c48400b80c8fd391f10cdd06a4f1d9c53ec490','family':'six-row trapezoidal SBP block; reflected outer block','nodes':n,'arithmetic':'exact rational/symbolic','sympy':s.__version__,'free_parameters':[str(z) for z in free], 'inner_rows':[[str(s.factor(D[i,j])) for j in range(m+1)] for i in range(m)],'quadratic_conditions_rows_1_through_6':[str(z) for z in conditions],'left_error_weighted_sum':'1/2','weighted_residual_identity':str(invariant),'boundary_moments':'INFEASIBLE','scalar_SBP':True,'Dirichlet_kernel':'UNRESOLVED: no moment-feasible member exists','zero_parameter_scalar_Dirichlet_rank':rank,'zero_parameter_old_kernel_detected':D0*old!=s.zeros(n,1),'full_STF_adjoint':'PASS by symbolic pairing identity in derivation','coercivity':'UNRESOLVED for a consistent member; none exists','physical_solves':0,'grid_sweeps':0}
Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
