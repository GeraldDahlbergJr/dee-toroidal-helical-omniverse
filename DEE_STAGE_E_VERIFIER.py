#!/usr/bin/env python3
"""Stage E exact local identities; no global existence/topology certification.
Run: python DEE_STAGE_E_VERIFIER.py (Python 3, SymPy).
Uses arbitrary two-coordinate functions in local isothermal meridional charts.
"""
import sympy as s

COUNT = 0

def check(label, value):
    global COUNT
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    for e in entries:
        residual = s.simplify(s.expand(e.doit()))
        if residual != 0:
            print('FAIL:', label, residual, flush=True)
            raise AssertionError(label)
    COUNT += 1
    print('PASS:', label, flush=True)


def connection(g, inv, xs):
    n = len(xs)
    return [[[s.simplify(sum(inv[i,k]*(g[k,j].diff(xs[l])+g[k,l].diff(xs[j])
                -g[j,l].diff(xs[k])) for k in range(n))/2)
              for l in range(n)] for j in range(n)] for i in range(n)]


def ricci(Ga, xs):
    n = len(xs)
    return s.Matrix(n,n,lambda i,j:s.simplify(
        sum(Ga[k][i][j].diff(xs[k])-Ga[k][i][k].diff(xs[j]) for k in range(n))
        +sum(Ga[k][k][l]*Ga[l][i][j]-Ga[k][j][l]*Ga[l][i][k]
             for k in range(n) for l in range(n))))


def main():
    print('SymPy', s.__version__, flush=True)
    t,r,z,p = xs = s.symbols('t r z phi', real=True)
    ys = (r,z)
    N,H,W,o,u,q1,q2 = [s.Function(n)(r,z) for n in
                       ('N','H','W','omega','rho','Theta','Psi')]
    a,b,c,lam,v = s.symbols('a b c lambda_V v', real=True)
    w=s.Matrix(s.symbols('w1 w2', real=True)); m=s.Matrix(s.symbols('m1 m2', real=True))
    q=s.Matrix([q1,q2]); U=w+o*m
    x=s.Symbol('x'); A=s.Matrix([[a,c],[c,b]])
    Mx=s.eye(2)+x*x*A; M=Mx.subs(x,u); Mp=Mx.diff(x).subs(x,u)
    Vx=lam*(x*x-v*v)**2/4; V=Vx.subs(x,u)
    g=s.Matrix([[-N*N+W*W*o*o,0,0,-W*W*o],
                [0,H*H,0,0],[0,0,H*H,0],[-W*W*o,0,0,W*W]])
    inv=s.Matrix([[-1/N**2,0,0,-o/N**2],[0,1/H**2,0,0],
                  [0,0,1/H**2,0],[-o/N**2,0,0,1/W**2-o**2/N**2]])
    check('inverse metric',g*inv-s.eye(4))
    check('determinant',g.det()+N*N*H**4*W*W)
    vol=N*H*H*W
    dQ=[w,q.diff(r),q.diff(z),m]
    C=s.Matrix(4,4,lambda i,j:u.diff(xs[i])*u.diff(xs[j])+(dQ[i].T*M*dQ[j])[0])
    P=(U.T*M*U)[0]/N**2
    X=sum(u.diff(y)**2+(q.diff(y).T*M*q.diff(y))[0] for y in ys)/H**2+(m.T*M*m)[0]/W**2
    L=(P-X)/2-V; T=C+g*L
    check('kinetic contraction and all temporal-angular cross terms',s.trace(inv*C)-X+P)
    # Independent inverse-metric variation, symmetric covariant kinetic C fixed.
    for i in range(4):
        for j in range(i,4):
            f=s.Symbol('metric_variation')
            variation=s.zeros(4); variation[i,j]=1; variation[j,i]=1
            gv=inv+f*variation
            lv=-s.trace(gv*C)/2-V
            varied=-2*s.diff(lv,f).subs(f,0)/(1 if i==j else 2)+g[i,j]*L
            check('metric variation T '+str((i,j)),varied-T[i,j])
    er=sum(s.diff(N*W*u.diff(y),y) for y in ys)/vol-Vx.diff(x).subs(x,u)
    er-=(sum((q.diff(y).T*Mp*q.diff(y))[0] for y in ys)/H**2
         +(m.T*Mp*m)[0]/W**2-(U.T*Mp*U)[0]/N**2)/2
    currents=[sum((inv[i,j]*M*dQ[j] for j in range(4)),s.zeros(2,1)) for i in range(4)]
    eq=sum((s.diff(vol*currents[i],xs[i]) for i in range(4)),s.zeros(2,1))/vol
    for f,e in [(u,er),(q1,eq[0]),(q2,eq[1])]:
        el=(s.diff(vol*L,f)-sum(s.diff(s.diff(vol*L,f.diff(y)),y) for y in ys))/vol
        check('action vs covariant equation '+str(f.func),el-e)
    check('temporal current',currents[0]+M*U/N**2)
    check('angular current',currents[3]-M*m/W**2+o*M*U/N**2)
    for k in range(2):
        check('meridional conserved current '+str(k),eq[k]-sum(s.diff(N*W*(M*q.diff(y))[k],y) for y in ys)/vol)
    Ga=connection(g,inv,xs); mixed=inv*T
    for j in range(4):
        div=sum(mixed[i,j].diff(xs[i]) for i in range(4))+sum(
            Ga[i][i][k]*mixed[k,j]-Ga[k][i][j]*mixed[i,k]
            for i in range(4) for k in range(4))
        check('off-shell stress divergence '+str(xs[j]),div-er*u.diff(xs[j])-(dQ[j].T*eq)[0])
    n=s.Matrix([1/N,0,0,o/N]); eps=(P+X)/2+V
    check('positive observer energy',(n.T*T*n)[0]-eps)
    jphi=-(U.T*M*m)[0]/N
    check('matter angular momentum',-(n.T*T)[3]-jphi)
    stress=sum(T[i,i]/([H**2,H**2,W**2][i-1]) for i in (1,2,3))
    check('stationary active mass',eps+stress-2*P+2*V)
    # Independent direct spacetime curvature and spatial ADM geometry.
    Ric=ricci(Ga,xs); R4=s.trace(inv*Ric)
    h=s.diag(H**2,H**2,W**2); hi=s.diag(1/H**2,1/H**2,1/W**2)
    spatial=(r,z,p); G3=connection(h,hi,spatial); R3ij=ricci(G3,spatial)
    K=s.zeros(3)
    for i,y in enumerate(ys): K[i,2]=K[2,i]=-W**2*o.diff(y)/(2*N)
    Ksq=s.trace(hi*K*hi*K)
    check('extrinsic curvature trace',s.trace(hi*K))
    check('extrinsic curvature norm',Ksq-W**2*sum(o.diff(y)**2 for y in ys)/(2*N**2*H**2))
    # Verify K from its defining stationary shift equation.
    beta=s.Matrix([0,0,-o]); bc=h*beta
    for i in range(3):
        for j in range(i,3):
            dk=(bc[j].diff(spatial[i])+bc[i].diff(spatial[j])
                -2*sum(G3[k][i][j]*bc[k] for k in range(3)))/(2*N)
            check('shift definition K '+str((i,j)),dk-K[i,j])
    Gnn=(n.T*Ric*n)[0]+R4/2
    check('Hamiltonian Gauss identity',Gnn-(s.trace(hi*R3ij)-Ksq)/2)
    Km=hi*K
    divK=[]
    for i in range(3):
        di=sum(Km[j,i].diff(spatial[j]) for j in range(3))+sum(
            G3[j][j][k]*Km[k,i]-G3[k][j][i]*Km[j,k] for j in range(3) for k in range(3))
        divK.append(di)
        check('momentum Codazzi sign '+str(i),di+(n.T*Ric)[i+1])
    drag=sum(s.diff(W**3*o.diff(y)/N,y) for y in ys)
    check('frame dragging divergence',divK[2]+drag/(2*H**2*W))
    lapN=sum(s.diff(W*N.diff(y),y) for y in ys)/(H**2*W)
    check('normal Ricci lapse identity',(n.T*Ric*n)[0]-lapN/N+Ksq)
    # Full spatial projected Ricci identity, including Lie_beta K.
    for i in range(3):
        for j in range(i,3):
            Hess=N.diff(spatial[i],spatial[j])-sum(G3[k][i][j]*N.diff(spatial[k]) for k in range(3))
            Lie=sum(beta[k]*K[i,j].diff(spatial[k])+K[k,j]*beta[k].diff(spatial[i])
                    +K[i,k]*beta[k].diff(spatial[j]) for k in range(3))
            expected=R3ij[i,j]-Hess/N-2*(K*hi*K)[i,j]+Lie/N
            check('spatial Einstein geometry '+str((i,j)),Ric[i+1,j+1]-expected)
    for i in (1,2):
        check('circular Ricci tA '+str(i),Ric[0,i])
        check('circular Ricci phiA '+str(i),Ric[3,i])
    flux=sum(s.diff(o*W**3*o.diff(y)/N,y) for y in ys)
    check('rotation Green identity',flux-o*drag-W**3*sum(o.diff(y)**2 for y in ys)/N)
    # Asymptotic algebra in the orthonormal frame, spatial gradients decayed.
    mv=Mx.subs(x,v); w2=(w.T*mv*w)[0]
    check('vacuum temporal energy coefficient',P.subs({u:v,N:1,o:0})/2-w2/2)
    check('vacuum temporal stress isotropy',L.subs({u:v,N:1,o:0,m[0]:0,m[1]:0,
        u.diff(r):0,u.diff(z):0,q1.diff(r):0,q1.diff(z):0,q2.diff(r):0,q2.diff(z):0})-w2/2)
    check('winding coefficient at zero amplitude',(m.T*Mx.subs(x,0)*m)[0]-(m.T*m)[0])
    check('kinetic determinant',Mx.det()-(1+(a+b)*x*x+(a*b-c*c)*x**4))
    # Direct flat null-wave stress: d_t Q=d_z Q=k; invariant zero, T not zero.
    k=s.Matrix(s.symbols('k1 k2',real=True))
    flat=s.diag(-1,1,1,1); nullgrad=[k,s.zeros(2,1),k,s.zeros(2,1)]
    cn=s.Matrix(4,4,lambda i,j:(nullgrad[i].T*M*nullgrad[j])[0])
    kn=s.trace(flat*cn); tn=cn-flat*kn/2
    check('null-wave kinetic cancellation',kn)
    check('null-wave energy retains positive quadratic',tn[0,0]-(k.T*M*k)[0])
    check('null-wave flux retains full stress',tn[0,2]-(k.T*M*k)[0])
    check('rotation source coercive when w=0',
          (U.T*M*m)[0].subs({w[0]:0,w[1]:0})-o*(m.T*M*m)[0])
    # Flat asymptotic angular momentum sign for omega=2 G J/R^3.
    R,G,J,theta=s.symbols('R G J theta',positive=True)
    omega_tail=2*G*J/R**3
    jsurface=-s.integrate((R*s.sin(theta))**2*s.diff(omega_tail,R)
                         *R**2*s.sin(theta),(theta,0,s.pi))*2*s.pi/(16*s.pi*G)
    check('ADM angular momentum tail normalization',jsurface-J)
    print(f'PASS: Stage E {COUNT} exact local checks. Global proofs are in the reports.',flush=True)

if __name__ == '__main__':
    main()
