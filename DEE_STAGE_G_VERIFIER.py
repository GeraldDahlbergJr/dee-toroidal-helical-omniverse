#!/usr/bin/env python3
"""Stage G exact local checks. Python 3 / SymPy, no numerical sampling.
General jets/densities and a documented noncircular curvature family.
Global ODE/continuation proofs and their assumptions are in the reports.
"""
import sys
import sympy as s

COUNT = 0


def check(name, expr):
    global COUNT
    entries = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    for entry in entries:
        result = s.factor(s.simplify(s.expand(entry.doit())))
        if result != 0:
            print('FAIL:', name, result, flush=True)
            raise AssertionError(name)
    COUNT += 1
    print('PASS:', name, flush=True)


def scalar_checks():
    N, f, om, V = s.symbols('N f Omega V', real=True, nonzero=True)
    # General h: verify block congruence before taking determinants.
    h = s.Matrix([[s.Symbol('h'+str(min(i,j))+str(max(i,j)))
                   for j in range(3)] for i in range(3)])
    beta = s.Matrix(s.symbols('b0:3', real=True))
    shear = s.eye(4)
    for i in range(3):
        shear[i+1, 0] = beta[i]
    diag = s.diag(-N**2, h)
    g = shear.T*diag*shear
    gi = s.zeros(4); gi[0,0] = -1/N**2
    hi = h.inv()
    for i in range(3):
        gi[0,i+1] = gi[i+1,0] = beta[i]/N**2
        for j in range(3):
            gi[i+1,j+1] = hi[i,j]-beta[i]*beta[j]/N**2
    check('general ADM inverse', g*gi-s.eye(4))
    check('ADM block congruence', g-s.BlockMatrix([
        [s.Matrix([-N**2+(beta.T*h*beta)[0]]), beta.T*h], [h*beta,h]]).as_explicit())
    check('unit triangular shear determinant', shear.det()-1)
    check('det g = -N squared det h via congruence', diag.det()+N**2*h.det())
    # Arbitrary scalar jets in an orthonormal spatial frame: tensorial projections.
    rr = s.Matrix(s.symbols('r0:3', real=True))
    dr = s.Matrix(s.symbols('x0:3', real=True))
    di = s.Matrix(s.symbols('y0:3', real=True))
    x,y = s.symbols('x y', real=True)
    zz=x+s.I*y; dz=dr+s.I*di
    bd=(beta.T*dz)[0]; q=-(beta.dot(rr))/N
    p=(s.I*om*zz-bd)/N
    norm=lambda a: s.expand(s.conjugate(a)*a)
    U=q*q+f*norm(p); X=rr.dot(rr)+f*(dr.dot(dr)+di.dot(di))
    lm=(U-X)/2-V
    gradR=s.Matrix([0,*rr]); gradX=s.Matrix([-om*y,*dr]); gradY=s.Matrix([om*x,*di])
    C=gradR*gradR.T+f*(gradX*gradX.T+gradY*gradY.T)
    gf=s.diag(-N**2,1,1,1)
    gf[0,0]+=beta.dot(beta)
    for i in range(3): gf[0,i+1]=gf[i+1,0]=beta[i]
    gif=gf.inv(); normal=s.Matrix([1,*(-beta)])/N
    check('Lorentz kinetic contraction',s.trace(gif*C)-(X-U))
    T=C+gf*lm
    check('coordinate Ttt',T[0,0]-(f*om**2*norm(zz)+gf[0,0]*lm))
    check('all coordinate Tti',s.Matrix([T[0,i+1]-f*om*s.im(s.conjugate(zz)*dz[i])-beta[i]*lm for i in range(3)]))
    check('positive ADM energy',(normal.T*T*normal)[0]-(U+X)/2-V)
    ji=s.Matrix([-sum(normal[mu]*T[mu,i+1] for mu in range(4)) for i in range(3)])
    jexpected=s.Matrix([-q*rr[i]-f*s.re(s.conjugate(p)*dz[i]) for i in range(3)])
    check('all ADM momentum projections',ji-jexpected)
    S=T[1:4,1:4]
    check('all spatial stresses',S-C[1:4,1:4]-s.eye(3)*lm)
    check('active mass epsilon plus S',(U+X)/2+V+s.trace(S)-2*(U-V))
    check('trace reversed spatial source',S-s.eye(3)*(s.trace(S)-((U+X)/2+V))/2-C[1:4,1:4]-s.eye(3)*V)
    check('normal phase kinetic shift signs', norm(p)-(om**2*norm(zz)+norm(bd)-2*om*s.im(s.conjugate(zz)*bd))/N**2)
    check('Komar rotation charge cancellation',N*U-beta.dot(ji)-f*om*s.im(s.conjugate(zz)*p))
    t=s.symbols('t',real=True)
    rot=s.Matrix([[s.cos(om*t),-s.sin(om*t)],[s.sin(om*t),s.cos(om*t)]])
    check('stationary stress common rotation',rot.T*rot-s.eye(2))
    # Metric variation of every symmetric component, unconstrained inverse jets.
    ij=[(i,j) for i in range(4) for j in range(i,4)]
    ss=s.symbols('g0:'+str(len(ij)))
    inv=s.zeros(4)
    for a,(i,j) in zip(ss,ij): inv[i,j]=inv[j,i]=a
    lag=-s.trace(inv*C)/2-V
    stress=s.Matrix(4,4,lambda i,j:C[i,j]+gf[i,j]*lag)
    for a,(i,j) in zip(ss,ij):
        check('metric variation stress '+str((i,j)),
              -2*s.diff(lag,a)/(1 if i==j else 2)+gf[i,j]*lag-stress[i,j])
    # General three-coordinate density differentiation, no suppressed shift jets.
    coords=s.symbols('X Y Z',real=True); t=s.symbols('t',real=True)
    fun=lambda name:s.Function(name)(*coords)
    nn,ff,vol=map(fun,('lapse','stiff','sqrt_h'))
    bv=[fun('beta'+str(i)) for i in range(3)]
    hh=s.Matrix(3,3,lambda i,j:fun('hi'+str(min(i,j))+str(max(i,j))))
    pp=hh-s.Matrix(bv)*s.Matrix(bv).T/nn**2
    weight=nn*vol; z=fun('zeta'); rho=fun('rho')
    field=s.exp(s.I*om*t)*z
    gti=s.Matrix(bv)/nn**2
    direct=s.diff(weight*ff*(-field.diff(t)/nn**2+sum(gti[i]*field.diff(coords[i]) for i in range(3))),t)
    direct+=sum(s.diff(weight*ff*(gti[i]*field.diff(t)+sum(pp[i,j]*field.diff(coords[j]) for j in range(3))),coords[i]) for i in range(3))
    reduced=sum(s.diff(weight*ff*sum(pp[i,j]*z.diff(coords[j]) for j in range(3)),coords[i]) for i in range(3))/weight
    reduced+=2*s.I*om*ff*sum(bv[i]*z.diff(coords[i]) for i in range(3))/nn**2
    reduced+=(s.I*om/weight*sum(s.diff(weight*ff*bv[i]/nn**2,coords[i]) for i in range(3))+om**2*ff/nn**2)*z
    check('general 3D covariant phase density reduction',direct/(weight*s.exp(s.I*om*t))-reduced)
    alternate=sum(s.diff(weight*ff*(sum(hh[i,j]*z.diff(coords[j]) for j in range(3))+bv[i]*(s.I*om*z-sum(bv[j]*z.diff(coords[j]) for j in range(3)))/nn**2),coords[i]) for i in range(3))
    alternate-=s.I*om*weight*ff*(s.I*om*z-sum(bv[j]*z.diff(coords[j]) for j in range(3)))/nn**2
    check('unexpanded phase current equivalence',alternate/weight-reduced)
    covrho=sum(s.diff(weight*sum(pp[i,j]*rho.diff(coords[j]) for j in range(3)),coords[i]) for i in range(3))/weight
    # Euler derivative of stationary amplitude kinetic term, using independent jets.
    rjets=s.symbols('rj0:3'); rvar=s.symbols('rvar'); a,lam,v=s.symbols('a lambda_V v',real=True)
    zre=[fun('zr'),fun('zi')]
    ztem=[-om*zre[1],om*zre[0]]
    zcontract=sum(sum(hh[i,j]*s.diff(u,coords[i])*s.diff(u,coords[j]) for i in range(3) for j in range(3))-(tau-sum(bv[i]*s.diff(u,coords[i]) for i in range(3)))**2/nn**2 for u,tau in zip(zre,ztem))
    kinetic=sum(pp[i,j]*rjets[i]*rjets[j] for i in range(3) for j in range(3))
    pot=lam*(rvar**2-v**2)**2/4
    density=weight*(-kinetic/2-(1+a*rvar**2)*zcontract/2-pot)
    euler=density.diff(rvar).subs(rvar,rho)-sum(s.diff(density.diff(rjets[i]).subs(dict(zip(rjets,[rho.diff(c) for c in coords]))),coords[i]) for i in range(3))
    check('amplitude action versus covariant reduction',euler/weight-covrho+lam*rho*(rho**2-v**2)+a*rho*zcontract)
    check('flat Helmholtz limit',reduced.subs({nn:1,ff:1,vol:1,**dict(zip(bv,[0]*3)),**{hh[i,j]:s.Integer(i==j) for i in range(3) for j in range(i,3)}}).doit()-sum(z.diff(c,2) for c in coords)-om**2*z)
    check('vacuum radial mass',s.diff(pot,rvar,2).subs(rvar,v)-2*lam*v**2)
    eta=s.symbols('eta')
    check('amplitude perturbation polynomial',s.diff(pot,rvar).subs(rvar,v+eta)-2*lam*v**2*eta-lam*(3*v*eta**2+eta**3))
    # Boundary product rules: arbitrary 1D representative of summed index identity.
    xx=coords[0]; weight=fun('weight'); A=fun('A'); F=fun('F'); B=fun('B')
    u=fun('u'); vv=fun('vv')
    def L(w):return s.diff(A*s.diff(w,xx),xx)/weight+s.I*om*(2*B*s.diff(w,xx)+s.diff(weight*B,xx)*w/weight)
    # vbar treated as independent; conjugating L reverses i.
    green=weight*(vv*L(u)-u*L(vv).subs(s.I,-s.I))
    flux=A*(vv*u.diff(xx)-u*vv.diff(xx))+2*s.I*om*weight*B*u*vv
    check('formal Hermitian Green boundary identity',green-s.diff(flux,xx))
    check('amplitude boundary product rule',s.diff(u*A*u.diff(xx),xx)-u*s.diff(A*u.diff(xx),xx)-A*u.diff(xx)**2)
    check('momentum integration product rule',s.diff(B*F,xx)-B*s.diff(F,xx)-F*s.diff(B,xx))
    # Real phase Green identity including the full linear frequency/shift term.
    zr,zi=fun('zr_boundary'),fun('zi_boundary')
    Cshift,Dfreq=fun('Cshift'),fun('Dfreq')
    # Functions here have real coefficients/values; explicitly conjugate through jets.
    op_re=s.diff(A*zr.diff(xx),xx)-om*(2*Cshift*zi.diff(xx)+Cshift.diff(xx)*zi)+om**2*Dfreq*zr
    op_im=s.diff(A*zi.diff(xx),xx)+om*(2*Cshift*zr.diff(xx)+Cshift.diff(xx)*zr)+om**2*Dfreq*zi
    flux=A*(zr*zr.diff(xx)+zi*zi.diff(xx))
    bulk=A*(zr.diff(xx)**2+zi.diff(xx)**2)+2*om*Cshift*(zr*zi.diff(xx)-zi*zr.diff(xx))-om**2*Dfreq*(zr**2+zi**2)
    check('phase real Green boundary and all shift terms',zr*op_re+zi*op_im-s.diff(flux,xx)+bulk)



def asymptotic_checks():
    e,m,k,om,fv,d=s.symbols('epsilon m k Omega f_v d',real=True, nonzero=True)
    N=1-m*e; H=1+m*e
    leading=lambda x:s.series(x,e,0,2).removeO()
    check('monopole volume density',leading(N*H**3)-(1+2*m*e))
    check('monopole spatial density cancellation',leading(N*H)-1)
    check('monopole frequency density',leading(H**3/N)-(1+4*m*e))
    check('isotropic optical frequency coefficient',leading(H**2/N**2)-(1+4*m*e))
    check('effective Coulomb potential sign',leading(om**2*(1-H**2/N**2))+4*m*om**2*e)
    R=s.symbols('R',positive=True)
    n=s.Function('N')(R); hh=s.Function('H')(R); ff=s.Function('f')(R)
    B=s.Function('B')(R); y=s.Function('y')(R); ell=s.symbols('ell',integer=True)
    # General conjugation independent of any asymptotic substitution.
    u=s.Function('u')(R); z=s.exp(-B/2)*u
    raw=z.diff(R,2)+(2/R+B.diff(R))*z.diff(R)+om**2*hh**2/n**2*z
    transformed=u.diff(R,2)+2*u.diff(R)/R+(om**2*hh**2/n**2-B.diff(R,2)/2-B.diff(R)/R-B.diff(R)**2/4)*u
    check('exact radial density conjugation',s.exp(B/2)*raw-transformed)
    zl=s.exp(-B/2)*y/R
    ode=zl.diff(R,2)+(2/R+B.diff(R))*zl.diff(R)+(om**2*hh**2/n**2-ell*(ell+1)/R**2)*zl
    expected=y.diff(R,2)+(om**2*hh**2/n**2-ell*(ell+1)/R**2-B.diff(R,2)/2-B.diff(R)/R-B.diff(R)**2/4)*y
    check('spherical harmonic radial equation',R*s.exp(B/2)*ode-expected)
    c=s.symbols('c',real=True)
    theta=k*R+c*s.log(R)/(2*k)
    for sign in (1,-1):
        w=s.exp(sign*s.I*theta)
        check('Coulomb comparison residual '+str(sign),(w.diff(R,2)+(k**2+c/R)*w)/w+c**2/(4*k**2*R**2)+sign*s.I*c/(2*k*R**2))
    wp,wm=s.exp(s.I*theta),s.exp(-s.I*theta)
    check('Coulomb comparison Wronskian',s.det(s.Matrix([[wp,wm],[wp.diff(R),wm.diff(R)]]))+2*s.I*theta.diff(R))
    eta=s.symbols('eta'); a,v=s.symbols('a v')
    check('asymptotic amplitude stiffness',(1+a*(v+eta)**2)-(1+a*v**2+2*a*v*eta+a*eta**2))


def curvature_checks():
    # Exact noncircular stationary family; all components, arbitrary functions x.
    t,x,y,z=coords=s.symbols('t x y z',real=True)
    N,H,W,b,c=[s.Function(n)(x) for n in ('N','H','W','b','c')]
    h=s.diag(1,H**2,W**2); hi=h.inv(); beta=s.Matrix([b,c,0])
    g=s.zeros(4);g[0,0]=-N**2+(beta.T*h*beta)[0]
    for i in range(3):
        g[0,i+1]=g[i+1,0]=(h*beta)[i]
        for j in range(3):g[i+1,j+1]=h[i,j]
    gi=g.inv()
    def connection(metric, inv, cc):
        dim=len(cc)
        return [[[s.factor(sum(inv[i,l]*(metric[l,j].diff(cc[k])+metric[l,k].diff(cc[j])-metric[j,k].diff(cc[l])) for l in range(dim))/2) for k in range(dim)] for j in range(dim)] for i in range(dim)]
    def ricci(ga,cc):
        dim=len(cc)
        return s.Matrix(dim,dim,lambda i,j:s.factor(sum(ga[k][i][j].diff(cc[k])-ga[k][i][k].diff(cc[j])+sum(ga[k][k][l]*ga[l][i][j]-ga[k][j][l]*ga[l][i][k] for l in range(dim)) for k in range(dim))))
    ga=connection(g,gi,coords); ric=ricci(ga,coords)
    cc=coords[1:]; gh=connection(h,hi,cc); rh=ricci(gh,cc)
    bc=h*beta
    K=s.Matrix(3,3,lambda i,j:s.factor((bc[j].diff(cc[i])+bc[i].diff(cc[j])-2*sum(gh[l][i][j]*bc[l] for l in range(3)))/(2*N)))
    tr=s.factor(s.trace(hi*K)); K2=s.factor(s.trace(hi*K*hi*K))
    dn=s.Matrix(3,3,lambda i,j:N.diff(cc[i],cc[j])-sum(gh[l][i][j]*N.diff(cc[l]) for l in range(3)))
    lie=s.Matrix(3,3,lambda i,j:sum(beta[k]*K[i,j].diff(cc[k])+K[k,j]*beta[k].diff(cc[i])+K[i,k]*beta[k].diff(cc[j]) for k in range(3)))
    normal=s.Matrix([1,-b,-c,0])/N
    Rnn=s.factor((normal.T*ric*normal)[0])
    R3=s.trace(hi*rh); R4=s.trace(gi*ric)
    check('nonmaximal K trace',tr-(b.diff(x)+b*(H.diff(x)/H+W.diff(x)/W))/N)
    check('Hamiltonian direct curvature Gauss',Rnn+R4/2-(R3+tr**2-K2)/2)
    A=hi*K-s.eye(3)*tr
    div=[]
    for i in range(3):
        div.append(sum(A[j,i].diff(cc[j])+sum(gh[j][j][l]*A[l,i]-gh[l][j][i]*A[j,l] for l in range(3)) for j in range(3)))
        check('momentum direct Codazzi '+str(i),sum(normal[mu]*ric[mu,i+1] for mu in range(4))+div[i])
    check('stationary normal Ricci lapse with beta DK',Rnn-s.trace(hi*dn)/N+K2+sum(beta[i]*tr.diff(cc[i]) for i in range(3))/N)
    spatial=rh+tr*K-2*K*hi*K-dn/N+lie/N
    for i in range(3):
        for j in range(i,3):check('spatial Ricci ADM '+str((i,j)),ric[i+1,j+1]-spatial[i,j])
    # Matter substitution in geometric constraints uses the verified projections.
    U,X,V,G=s.symbols('U X V G')
    eps=(U+X)/2+V; S=(3*U-X)/2-3*V
    check('lapse source coefficient',4*s.pi*G*(eps+S)-8*s.pi*G*(U-V))
    check('Hamiltonian matter coefficient',16*s.pi*G*eps-8*s.pi*G*(U+X+2*V))


def main():
    print('SymPy',s.__version__,flush=True)
    scalar_checks()
    asymptotic_checks()
    curvature_checks()
    print(f'PASS: Stage G {COUNT} exact local checks; global proofs are conditional as stated in the reports.',flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('FAIL: Stage G execution:',type(exc).__name__,str(exc),flush=True)
        sys.exit(1)
