#!/usr/bin/env python3
"""Exact identities only; global hypotheses and proofs are in Stage D reports."""
import sympy as s

COUNT = 0

def check(label, expr):
    global COUNT
    es = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    assert all(s.simplify(s.expand(e.doit())) == 0 for e in es), label
    COUNT += 1
    print('PASS:', label)


def main():
    print('SymPy', s.__version__)
    t, r, p, z = xs = s.symbols('t r phi z', real=True)
    N, H, W, u, q1, q2 = [s.Function(n)(r, z) for n in ('N','H','W','rho','Theta','Psi')]
    q = s.Matrix([q1,q2]); x = s.Symbol('x')
    a,b,c,lam,v = s.symbols('a b c lambda_V v', real=True)
    Mx = s.eye(2)+x*x*s.Matrix([[a,c],[c,b]])
    Vx = lam*(x*x-v*v)**2/4
    M,V = Mx.subs(x,u),Vx.subs(x,u)
    F=N*W; S=F*H**2
    kinetic=sum(u.diff(y)**2+(q.diff(y).T*M*q.diff(y))[0] for y in (r,z))/H**2
    lag=-kinetic/2-V
    er=sum(s.diff(F*u.diff(y),y) for y in (r,z))/S-Vx.diff(x).subs(x,u)
    er-=sum((q.diff(y).T*Mx.diff(x).subs(x,u)*q.diff(y))[0] for y in (r,z))/(2*H**2)
    eq=sum((s.diff(F*M*q.diff(y),y) for y in (r,z)),s.zeros(2,1))/S
    density=S*lag
    for f,e in [(u,er),(q1,eq[0]),(q2,eq[1])]:
        el=(s.diff(density,f)-sum(s.diff(s.diff(density,f.diff(y)),y) for y in (r,z)))/S
        check('variational vs covariant '+str(f.func),el-e)
    g=s.diag(-N**2,H**2,W**2,H**2); inv=s.diag(-1/N**2,1/H**2,1/W**2,1/H**2)
    C=s.Matrix(4,4,lambda i,j:u.diff(xs[i])*u.diff(xs[j])+(q.diff(xs[i]).T*M*q.diff(xs[j]))[0])
    T=C+g*lag
    Ga=[[[sum(inv[i,k]*(g[k,j].diff(xs[l])+g[k,l].diff(xs[j])-g[j,l].diff(xs[k])) for k in range(4))/2 for l in range(4)] for j in range(4)] for i in range(4)]
    mixed=inv*T
    for j in range(4):
        div=sum(mixed[i,j].diff(xs[i]) for i in range(4))+sum(Ga[i][i][k]*mixed[k,j]-Ga[k][i][j]*mixed[i,k] for i in range(4) for k in range(4))
        check('off-shell stress divergence '+str(xs[j]),div-er*u.diff(xs[j])-(q.diff(xs[j]).T*eq)[0])
    Ric=s.Matrix(4,4,lambda i,j:sum(Ga[k][i][j].diff(xs[k])-Ga[k][i][k].diff(xs[j]) for k in range(4))+sum(Ga[k][k][l]*Ga[l][i][j]-Ga[k][j][l]*Ga[l][i][k] for k in range(4) for l in range(4)))
    scalar=s.trace(inv*Ric)
    B=s.log(H)
    R3=-2/H**2*(B.diff(r,2)+B.diff(z,2)+(W.diff(r,2)+W.diff(z,2))/W)
    lapN=sum(s.diff(W*N.diff(y),y) for y in (r,z))/(H**2*W)
    check('Hamiltonian geometry G_nn=R3/2',(Ric[0,0]+N**2*scalar/2)/N**2-R3/2)
    check('lapse geometry R_tt=N D2N',Ric[0,0]-N*lapN)
    for i in (1,2,3):
        check('momentum geometry '+str(i),Ric[0,i])
    hess=N.diff(r,z)-B.diff(r)*N.diff(z)-B.diff(z)*N.diff(r)
    r3rz=-W.diff(r,z)/W+(B.diff(z)*W.diff(r)+B.diff(r)*W.diff(z))/W
    check('mixed Einstein geometry',Ric[1,3]-r3rz+hess/N)
    eps=kinetic/2+V; stress=sum(T[i,i]*inv[i,i] for i in (1,2,3))
    check('static active mass eps+traceS=-2V',eps+stress+2*V)
    check('energy',T[0,0]/N**2-eps)
    # Product-rule identities including both phase currents, with no on-shell substitution.
    boundary=sum(s.diff(F*(q.T*M*q.diff(y))[0],y) for y in (r,z))
    phasebulk=sum((q.diff(y).T*M*q.diff(y))[0] for y in (r,z))*F
    check('phase integral identity',boundary-phasebulk-S*(q.T*eq)[0])
    ampflux=sum(s.diff(F*(u-v)*u.diff(y),y) for y in (r,z))
    check('amplitude integration identity',ampflux-F*sum(u.diff(y)**2 for y in (r,z))-(u-v)*sum(s.diff(F*u.diff(y),y) for y in (r,z)))
    check('positive-domain potential factor',(x-v)*Vx.diff(x)-lam*x*(x+v)*(x-v)**2)
    for f in (u,q1,q2):
        check('mixed partial closure '+str(f.func),f.diff(r,z)-f.diff(z,r))
    k=s.Function('k')(r,z)
    check('axial gradient cannot be specified alone',s.diff(k*z,r)-z*k.diff(r))
    for sign in (1,-1):
        check('vacuum V and Vprime '+str(sign),s.Matrix([Vx.subs(x,sign*v),Vx.diff(x).subs(x,sign*v)]))
    check('nonvacuum zero amplitude',Vx.subs(x,0)-lam*v**4/4)
    check('kinetic determinant',Mx.det()-(1+(a+b)*x*x+(a*b-c*c)*x**4))
    check('vacuum amplitude mass',Vx.diff(x,2).subs(x,v)-2*lam*v*v)
    print(f'PASS: Stage D {COUNT} exact checks. Global hypotheses are not symbolic conclusions.')

if __name__ == '__main__':
    main()
