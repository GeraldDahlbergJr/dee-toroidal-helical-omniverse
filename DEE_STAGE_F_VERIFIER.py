#!/usr/bin/env python3
"""Exact local Stage F checks. Python 3 + SymPy; no numerical sampling.
Completeness, topology and global energy theorems are analytical, in the reports.
One check may contain multiple component identities; totals count named checks.
"""
import sympy as s

COUNT = 0

def check(label, expr):
    global COUNT
    entries = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    for entry in entries:
        residual = s.factor(s.simplify(entry.doit()))
        if residual != 0:
            print('FAIL:', label, residual, flush=True)
            raise AssertionError(label)
    COUNT += 1
    print('PASS:', label, flush=True)


def main():
    print('SymPy', s.__version__, flush=True)
    r, x, y = q = s.symbols('rho x y', real=True)
    a, b, lam, v = s.symbols('alpha beta lambda_V v', real=True)
    f, h = 1+a*r*r, 1+b*r*r
    g = s.diag(1, f, h); gi = g.inv()
    ga = [[[s.simplify(sum(gi[i,l]*(g[l,j].diff(q[k])+g[l,k].diff(q[j])
            -g[j,k].diff(q[l])) for l in range(3))/2)
            for k in range(3)] for j in range(3)] for i in range(3)]
    expected = {(0,1,1):-a*r, (0,2,2):-b*r,
                (1,0,1):a*r/f, (1,1,0):a*r/f,
                (2,0,2):b*r/h, (2,2,0):b*r/h}
    check('all target Christoffel components', s.Matrix([
        ga[i][j][k]-expected.get((i,j,k),0)
        for i in range(3) for j in range(3) for k in range(3)]))
    ric = s.Matrix(3,3,lambda i,j:sum(ga[k][i][j].diff(q[k])-ga[k][i][k].diff(q[j])
        +sum(ga[k][k][l]*ga[l][i][j]-ga[k][j][l]*ga[l][i][k]
             for l in range(3)) for k in range(3)))
    scalar = s.factor(s.trace(gi*ric))
    formula = -2*(a/f**2+b/h**2+a*b*r*r/(f*h))
    check('direct Ricci scalar', scalar-formula)
    numerator = a+b+5*a*b*r*r+2*a*b*(a+b)*r**4+a*a*b*b*r**6
    check('curvature rational numerator', scalar*f*f*h*h+2*numerator)
    check('rank one curvature', scalar.subs(b,0)+2*a/f**2)
    check('opposite eigenvalues are not flat', numerator.subs(b,-a)+5*a*a*r*r-a**4*r**6)
    check('repeated curvature', scalar.subs(b,a)+2*a*(2+a*r*r)/f**2)
    check('only flat parameter candidate from numerator coefficients',
          s.Matrix([numerator.subs(r,0)-(a+b), s.expand(numerator).coeff(r,2)-5*a*b]))
    def lie(X, metric=g):
        return s.Matrix(3,3,lambda i,j:sum(X[k]*metric[i,j].diff(q[k])
            +metric[k,j]*X[k].diff(q[i])+metric[i,k]*X[k].diff(q[j]) for k in range(3)))
    U, P, H = [s.Function(n)(*q) for n in ('U','P','H')]
    X = s.Matrix([U,P,H]); lk = lie(X)
    ke = s.Matrix([[2*U.diff(r), U.diff(x)+f*P.diff(r), U.diff(y)+h*H.diff(r)],
         [U.diff(x)+f*P.diff(r),2*a*r*U+2*f*P.diff(x), f*P.diff(y)+h*H.diff(x)],
         [U.diff(y)+h*H.diff(r),f*P.diff(y)+h*H.diff(x),2*b*r*U+2*h*H.diff(y)]])
    check('six unreduced Killing equations', lk-ke)
    ex, ey = s.Matrix([0,1,0]), s.Matrix([0,0,1])
    R = s.Matrix([0,-y,x]); er = s.Matrix([1,0,0])
    check('translation x', lie(ex)); check('translation y', lie(ey))
    check('rotation obstruction for distinct eigenvalues',lie(R)-s.Matrix([[0,0,0],[0,0,h-f],[0,h-f,0]]))
    check('repeated eigenvalue rotation',lie(R).subs(b,a))
    c, d, om = s.symbols('c d Omega', real=True)
    gen = s.Matrix([0,c-om*y,d+om*x])
    check('general affine phase generator obstruction',lie(gen)[1,2]-om*(b-a)*r*r)
    # Flat E(3): all generators, then every pairwise bracket and closure.
    flatbasis = [er, ex, ey]
    names = ['T_rho','T_x','T_y']
    for i,j in ((0,1),(0,2),(1,2)):
        z=s.zeros(3,1); z[j]=q[i]; z[i]=-q[j]
        flatbasis.append(z); names.append('L'+str(i)+str(j))
    for name,z in zip(names,flatbasis):
        check('flat Killing '+name,lie(z,s.eye(3)))
    def bracket(z,w):
        return w.jacobian(q)*z-z.jacobian(q)*w
    for i in range(6):
        for j in range(i+1,6):
            z=bracket(flatbasis[i],flatbasis[j]); B=z.jacobian(q)
            const=z.subs({r:0,x:0,y:0})
            reconstruction=const+sum((B[l,k]*flatbasis[3+idx]
                 for idx,(k,l) in enumerate(((0,1),(0,2),(1,2)))),s.zeros(3,1))
            check('flat algebra closure '+names[i]+','+names[j],
                  s.Matrix(list(B+B.T)+list(z-reconstruction)))
    check('phase translations commute',bracket(ex,ey))
    check('[R,T_x]=-T_y',bracket(R,ex)+ey)
    check('[R,T_y]=T_x',bracket(R,ey)-ex)
    V=lam*(r*r-v*v)**2/4
    for name,z in zip(names,flatbasis):
        check('potential action '+name, (s.Matrix([V.diff(zq) for zq in q]).dot(z))-V.diff(r)*z[0])
    check('phase affine potential preservation',V.diff(r)*gen[0])
    check('mixed flat rotation potential breaking', V.diff(r)*flatbasis[3][0]+lam*r*x*(r*r-v*v))
    # General jets in a Lorentz orthonormal spacetime frame. Tensorial identities
    # extend covariantly by spacetime normal coordinates, not by assuming flat g.
    signs=(-1,1,1,1)
    D=[s.Matrix(s.symbols('d'+str(mu)+'_0:3')) for mu in range(4)]
    B=[s.Matrix(s.symbols('b'+str(mu)+'_0:3')) for mu in range(4)]
    box=sum((signs[mu]*B[mu] for mu in range(4)),s.zeros(3,1))
    E=g*box+s.Matrix([sum(signs[mu]*g[i,l]*ga[l][j][k]*D[mu][j]*D[mu][k]
             for mu in range(4) for l in range(3) for j in range(3) for k in range(3))
             -V.diff(q[i]) for i in range(3)])
    def noether_residual(z):
        cov=g*z
        div=sum(signs[mu]*(cov.dot(B[mu])+(D[mu].T*cov.jacobian(q)*D[mu])[0]) for mu in range(4))
        return div-z.dot(E)-z[0]*V.diff(r)-sum(signs[mu]*(D[mu].T*lie(z)*D[mu])[0]/2 for mu in range(4))
    check('off-shell arbitrary vector Noether identity',noether_residual(X))
    for name,z in [('T_x',ex),('T_y',ey),('R',R),('T_rho',er),('L01',flatbasis[3])]:
        check('current contraction and sign '+name,noether_residual(z))
    check('rotation current components', (g*R)-s.Matrix([0,-f*y,h*x]))
    # Time-equivariant pullback: differentiate G and both scalar derivative legs.
    Lg=lie(X); DX=X.jacobian(q)
    C=s.Matrix(4,4,lambda mu,nu:(D[mu].T*g*D[nu])[0])
    dC=s.Matrix(4,4,lambda mu,nu:sum(X[i]*C[mu,nu].diff(q[i]) for i in range(3))
         +((DX*D[mu]).T*g*D[nu])[0]+(D[mu].T*g*DX*D[nu])[0])
    pull=s.Matrix(4,4,lambda mu,nu:(D[mu].T*Lg*D[nu])[0])
    check('stationary stress pullback kinetic identity',dC-pull)
    eta=s.diag(*signs)
    check('full stress Lie identity', dC-eta*(s.trace(eta*dC)/2+X[0]*V.diff(r))
          -(pull-eta*(s.trace(eta*pull)/2+X[0]*V.diff(r))))
    for name,z,subs in [('translations',gen,{om:0}),('isotropic rotation',R,{b:a})]:
        check('stationary full action generator '+name,lie(z).subs(subs))
    kinetic=sum(signs[mu]*C[mu,mu] for mu in range(4))
    T=C+eta*(-kinetic/2-V)
    eps=(D[0].T*g*D[0])[0]/2+sum((D[mu].T*g*D[mu])[0]/2 for mu in (1,2,3))+V
    check('positive observer energy',T[0,0]-eps)
    N=s.symbols('N',positive=True); beta=s.symbols('shift1:4')
    nd=(gen-sum((beta[i]*D[i+1] for i in range(3)),s.zeros(3,1)))/N
    check('ADM normal kinetic expansion', (nd.T*g*nd)[0]-(
        (gen.T*g*gen)[0]-2*sum(beta[i]*(gen.T*g*D[i+1])[0] for i in range(3))
        +sum(beta[i]*beta[j]*(D[i+1].T*g*D[j+1])[0] for i in range(3) for j in range(3)))/N**2)
    xc,yc=s.symbols('x_c y_c',real=True)
    centered=s.Matrix([0,-om*(y-yc),om*(x-xc)])
    check('centered rotation norm',(centered.T*g*centered)[0].subs(b,a)-om**2*f*((x-xc)**2+(y-yc)**2))
    check('fixed vacuum centered rotation',centered.subs({r:v,x:xc,y:yc}))
    check('translation vacuum norm',(gen.T*g*gen)[0].subs({om:0,r:v})-((1+a*v*v)*c*c+(1+b*v*v)*d*d))
    check('rotation center solves affine fixed point',gen.subs({x:-d/om,y:c/om}))
    check('Stage E additive normal velocity',(gen/N).subs(om,0)-s.Matrix([0,c/N,d/N]))
    check('Stage E M(0) stiffness',g.subs(r,0)-s.eye(3))
    check('vacuum radial mass',V.diff(r,2).subs(r,v)-2*lam*v*v)
    check('vacuum phase masses',s.Matrix([V.diff(x,2),V.diff(y,2)]))
    # Explicit equivariant flat-space wave equation and smooth kinematic profile.
    t,z,Omega,L=s.symbols('t z Omega L',real=True)
    amp=s.Function('u')(x,y,z)
    Z=s.exp(s.I*Omega*t)*amp
    check('massless rotating asymptotic Helmholtz operator',
          -Z.diff(t,2)+sum(Z.diff(w,2) for w in (x,y,z))
          -s.exp(s.I*Omega*t)*(sum(amp.diff(w,2) for w in (x,y,z))+Omega**2*amp))
    profile=s.exp(-(x*x+y*y+z*z))
    check('Gaussian profile is not a free harmonic solution',
          sum(profile.diff(w,2) for w in (x,y,z))+Omega**2*profile
          -(4*(x*x+y*y+z*z)-6+Omega**2)*profile)
    # Reduced stationary phase equation with arbitrary stationary coefficients.
    Hvol=s.Function('Hvol')(x,y); stiff=s.Function('stiff')(x,y)
    coords=(x,y)
    tt=s.Function('gtt')(x,y)
    ti=[s.Function('gt'+str(i))(x,y) for i in range(2)]
    spatial=s.Matrix([[s.Function('gxx')(x,y),s.Function('gxy')(x,y)],
                      [s.Function('gxy')(x,y),s.Function('gyy')(x,y)]])
    zz=s.Function('zeta')(x,y); phase=s.exp(s.I*Omega*t); field=phase*zz
    direct=s.diff(Hvol*stiff*(tt*field.diff(t)+sum(ti[i]*field.diff(coords[i]) for i in range(2))),t)
    direct+=sum(s.diff(Hvol*stiff*(ti[i]*field.diff(t)+sum(spatial[i,j]*field.diff(coords[j])
                    for j in range(2))),coords[i]) for i in range(2))
    reduced=sum(s.diff(Hvol*stiff*sum(spatial[i,j]*zz.diff(coords[j]) for j in range(2)),coords[i])
                for i in range(2))/Hvol
    reduced+=2*s.I*Omega*stiff*sum(ti[i]*zz.diff(coords[i]) for i in range(2))
    reduced+=(s.I*Omega*sum(s.diff(Hvol*stiff*ti[i],coords[i]) for i in range(2))/Hvol
               -Omega**2*stiff*tt)*zz
    check('noncircular stationary reduced phase equation',direct/(Hvol*phase)-reduced)
    # Subsidiary circular ansatz uses the same two original Cartesian phases.
    phi=s.symbols('phi',real=True); n=s.symbols('n',integer=True)
    drag,W=s.symbols('omega W',real=True); F=s.Function('F')(x,y)
    doublet=s.Matrix([F*s.cos(Omega*t+n*phi),F*s.sin(Omega*t+n*phi)])
    rot=s.Matrix([[0,-1],[1,0]])
    pi=(doublet.diff(t)+drag*doublet.diff(phi))/N
    check('circular residual normal derivative',pi-(Omega+drag*n)*rot*doublet/N)
    check('circular residual angular momentum',-stiff*pi.dot(doublet.diff(phi))
          +stiff*(Omega+drag*n)*n*F**2/N)
    check('circular residual temporal energy',stiff*pi.dot(pi)/2
          -stiff*(Omega+drag*n)**2*F**2/(2*N**2))
    check('circular residual angular energy',stiff*doublet.diff(phi).dot(doublet.diff(phi))/(2*W**2)
          -stiff*n**2*F**2/(2*W**2))
    check('circular residual mixed compatibility',s.Matrix([
          doublet.diff(t).dot(doublet.diff(x)),doublet.diff(phi).dot(doublet.diff(x))]))
    # Exact Gaussian integrability; this does not enforce field equations.
    gaussian_norm=s.integrate(s.exp(-2*x*x),(x,-s.oo,s.oo))**3
    check('finite Gaussian spatial L2 norm',gaussian_norm-(s.pi/2)**s.Rational(3,2))
    # Target lattice descent obstruction is algebraic; topology is analytical.
    ell1,ell2=s.symbols('ell1 ell2',real=True)
    check('rotation lattice discrepancy',gen.subs({x:x+ell1,y:y+ell2},simultaneous=True)
          -gen-s.Matrix([0,-om*ell2,om*ell1]))
    print(f'PASS: Stage F {COUNT} exact local checks; global proofs and limitations are in the reports.',flush=True)

if __name__ == '__main__':
    main()
