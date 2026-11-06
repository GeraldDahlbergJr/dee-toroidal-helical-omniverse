"""Symbolically verify continuum Hessian and longitudinal-divergence identities.
No finite differences, evaluator imports, baseline arrays, or solves.
"""
import json
from pathlib import Path
import sympy as s
r,t,p=s.symbols('r theta phi',positive=True,real=True)
f=s.Function('f')(r,t,p);q=1+r*s.cos(t)
er=s.Matrix([s.cos(t)*s.cos(p),s.cos(t)*s.sin(p),s.sin(t)])
et=s.Matrix([-s.sin(t)*s.cos(p),-s.sin(t)*s.sin(p),s.cos(t)])
ep=s.Matrix([-s.sin(p),s.cos(p),0])
def D(v,i):return er[i]*s.diff(v,r)+et[i]/r*s.diff(v,t)+ep[i]/q*s.diff(v,p)
g=s.Matrix([D(f,i) for i in range(3)])
hrr=s.diff(f,r,2);htt=s.diff(f,t,2)/r**2+s.diff(f,r)/r
hpp=s.diff(f,p,2)/q**2+s.cos(t)*s.diff(f,r)/q-s.sin(t)*s.diff(f,t)/(r*q)
hrt=s.diff(f,r,t)/r-s.diff(f,t)/r**2
hrp=s.diff(f,r,p)/q-s.cos(t)*s.diff(f,p)/q**2
htp=s.diff(f,t,p)/(r*q)+s.sin(t)*s.diff(f,p)/q**2
H=er*er.T*hrr+et*et.T*htt+ep*ep.T*hpp+(er*et.T+et*er.T)*hrt+(er*ep.T+ep*er.T)*hrp+(et*ep.T+ep*et.T)*htp
checks=[]
# The operators applied to each Cartesian gradient component include frame derivatives.
for i in range(3):
    for j in range(3):
        diff=s.expand(D(g[i],j)-H[i,j])
        coeffs=s.Poly(diff,*sorted(diff.atoms(s.Derivative),key=str)).coeffs() if diff else [0]
        ok=all(s.trigsimp(s.factor(v))==0 for v in coeffs)
        assert ok,(i,j)
        checks.append({'identity':f'H_cartesian_{i}{j}','passed':ok})
x=s.symbols('x y z');w=s.Matrix([s.Function('W'+str(i))(*x) for i in range(3)])
jac=w.jacobian(x);L=jac+jac.T-s.eye(3)*s.trace(jac)*s.Rational(2,3)
for i in range(3):
    lhs=sum(s.diff(L[i,j],x[j]) for j in range(3))
    rhs=sum(s.diff(w[i],c,2) for c in x)+s.diff(s.trace(jac),x[i])/3
    ok=s.simplify(lhs-rhs)==0;assert ok
    checks.append({'identity':f'div_LW_{i}=lap_Wi+grad_div_W_i/3','passed':ok})
out={'status':'PASS','method':'Exact symbolic chain rule for arbitrary smooth f and W; no discrete derivative or solver equation imports','checks':checks}
Path(__file__).with_name('continuum_identity_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
