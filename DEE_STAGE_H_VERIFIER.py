#!/usr/bin/env python3
"""Stage H: bounded exact local algebra; no global spectral theorem is tested.
Run from any directory with Python 3 and SymPy. No files are written.
Checks count named identities, not matrix entries or analytical assertions.
"""
import sys
import sympy as s

COUNT = 0

def check(name, expression):
    global COUNT
    entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    for entry in entries:
        answer = s.factor(s.expand(entry.doit()))
        if answer != 0:
            answer = s.simplify(answer)
        if answer != 0:
            raise AssertionError(f'{name}: {answer}')
    COUNT += 1
    print('PASS:', name, flush=True)


def main():
    print('SymPy', s.__version__, flush=True)
    x = s.symbols('x0:3', real=True)
    fun = lambda n: s.Function(n)(*x)
    div = lambda q: sum(s.diff(q[i], x[i]) for i in range(3))
    grad = lambda q: s.Matrix([s.diff(q, xx) for xx in x])
    lap = lambda q: div(grad(q))
    N, f, om = s.symbols('N f Omega', real=True, nonzero=True)
    h = s.Matrix(3, 3, lambda i,j: s.Symbol(f'h{min(i,j)}{max(i,j)}'))
    b = s.Matrix(s.symbols('b0:3', real=True))
    hi = h.inv()
    shear = s.eye(4)
    shear[1:4,0] = b
    g = shear.T*s.diag(-N**2,h)*shear
    gi = s.zeros(4)
    gi[0,0] = -1/N**2
    gi[0,1:4] = b.T/N**2
    gi[1:4,0] = b/N**2
    gi[1:4,1:4] = hi-b*b.T/N**2
    check('general ADM inverse', g*gi-s.eye(4))
    check('ADM determinant by unit shear', shear.det()-1)
    check('ADM diagonal-block determinant', s.diag(-N**2,h).det()+N**2*h.det())
    xi = s.Matrix(s.symbols('xi0:3', real=True))
    check('principal quadratic form', (xi.T*(hi-b*b.T/N**2)*xi)[0]-(xi.T*hi*xi)[0]+(b.dot(xi))**2/N**2)
    lam = s.symbols('lambda')
    check('orthonormal principal eigenvalues', ((s.eye(3)-b*b.T/N**2)-lam*s.eye(3)).det()-(1-lam)**2*(1-lam-b.dot(b)/N**2))

    nn, ff, w, z = map(fun, ('N','f','w','z'))
    bv = s.Matrix([fun('beta'+str(i)) for i in range(3)])
    hh = s.Matrix(3,3,lambda i,j: fun(f'hi{min(i,j)}{max(i,j)}'))
    pp = hh-bv*bv.T/nn**2
    t = s.symbols('t', real=True)
    Z = s.exp(s.I*om*t)*z
    direct = s.diff(w*ff*(-s.diff(Z,t)+bv.dot(grad(Z)))/nn**2,t)
    direct += div(w*ff*(pp*grad(Z)+bv*s.diff(Z,t)/nn**2))
    L = div(w*ff*pp*grad(z))/w + 2*s.I*om*ff*bv.dot(grad(z))/nn**2
    L += (s.I*om*div(w*ff*bv/nn**2)/w+om**2*ff/nn**2)*z
    check('general stationary covariant phase reduction', direct/s.exp(s.I*om*t)/w-L)
    check('density pencil identity', w*L-div(w*ff*pp*grad(z))-s.I*om*(2*w*ff*bv.dot(grad(z))/nn**2+div(w*ff*bv/nn**2)*z)-om**2*w*ff*z/nn**2)
    # u and vbar are independent functions; conjugation reverses i only.
    u, vb = fun('u'), fun('vbar')
    aa = s.Matrix(3,3,lambda i,j: fun(f'A{min(i,j)}{max(i,j)}'))
    bb = s.Matrix([fun('B'+str(i)) for i in range(3)])
    cc = fun('C')
    def pencil(q, sign=1):
        return div(aa*grad(q))+sign*s.I*om*(2*bb.dot(grad(q))+div(bb)*q)+om**2*cc*q
    flux = vb*aa*grad(u)-u*aa*grad(vb)+2*s.I*om*bb*u*vb
    check('full anisotropic Green flux', vb*pencil(u)-u*pencil(vb,-1)-div(flux))
    check('shift skew transport product rule', vb*(2*bb.dot(grad(u))+div(bb)*u)+u*(2*bb.dot(grad(vb))+div(bb)*vb)-2*div(bb*u*vb))

    # First variation in epsilon: angular dependence is allowed in all jets.
    e = s.symbols('epsilon')
    n, F = fun('n'), fun('F')
    H = s.Matrix(3,3,lambda i,j:fun(f'H{min(i,j)}{max(i,j)}'))
    T = s.trace(H)
    d = n+T/2
    A1 = (d+F)*s.eye(3)-H
    C1 = F+T/2-n
    check('first variation inverse spatial metric', (s.eye(3)+e*H)*(s.eye(3)-e*H)-s.eye(3)-e**2*(-H*H))
    check('first variation spatial determinant', s.diff((s.eye(3)+e*H).det(),e).subs(e,0)-T)
    check('first variation volume sqrt determinant', s.diff(s.sqrt((s.eye(3)+e*H).det()),e).subs(e,0)-T/2)
    check('first variation spacetime density', s.diff((1+e*n)*(1+e*T/2),e).subs(e,0)-d)
    check('first variation spatial density', s.diff((1+e*d)*(1+e*F)*(s.eye(3)-e*H),e).subs(e,0)-A1)
    check('first variation frequency density', s.diff((1+e*d)*(1+e*F)/(1+e*n)**2,e).subs(e,0)-C1)
    check('first variation shift density', s.diff((1+e*d)*(1+e*F)*e*bv/(1+e*n)**2,e).subs(e,0)-bv)
    check('first variation shift squared principal term', s.diff(e**2*bv*bv.T/(1+e*n)**2,e).subs(e,0))
    density1 = div(A1*grad(z))+om**2*C1*z
    divided1 = -(H* s.hessian(z,x)).trace() + (grad(n+T/2)-s.Matrix([sum(s.diff(H[i,j],x[i]) for i in range(3)) for j in range(3)])).dot(grad(z))-2*om**2*n*z
    check('metric differential versus scalar correction', density1-(d+F)*(lap(z)+om**2*z)-divided1-grad(F).dot(grad(z)))
    rho, v, a, eta, fv = s.symbols('rho v a eta f_v', real=True)
    check('amplitude stiffness exact expansion', 1+a*(v+eta)**2-(1+a*v**2+2*a*v*eta+a*eta**2))
    ffj = fun('stiffness'); Njj=fun('lapse'); vol=fun('sqrt_h')
    check('exact division removes multiplicative stiffness', div(vol*Njj*ffj*pp*grad(z))/(vol*Njj*ffj)-div(vol*Njj*pp*grad(z))/(vol*Njj)-grad(ffj).dot(pp*grad(z))/ffj)
    m, R = s.symbols('m R', real=True, positive=True)
    check('isotropic mass spatial density cancellation', ((-m+3*m)*s.eye(3)-2*m*s.eye(3)))
    check('isotropic mass Coulomb frequency', (3*m-(-m))-4*m)
    # Areal Cartesian Schwarzschild first jet H_ij=2m x_i x_j/r^3.
    r = s.sqrt(sum(xx**2 for xx in x)); nh=s.Matrix(x)/r
    Har = 2*m*nh*nh.T/r
    check('areal gauge trace', s.trace(Har)-2*m/r)
    check('areal gauge radial principal correction', (nh.T*(((-m/r+s.trace(Har)/2)*s.eye(3))-Har)*nh)[0]+2*m/r)
    gauge = m*nh
    check('areal versus isotropic spatial gauge', 2*m*s.eye(3)/r-Har-s.Matrix(3,3,lambda i,j:s.diff(gauge[i],x[j])+s.diff(gauge[j],x[i])))
    # Linearized stationary constraints: tensor gauge cancellation.
    Ric = lambda q:s.Matrix(3,3,lambda i,j:(sum(s.diff(q[k,j],x[k],x[i])+s.diff(q[k,i],x[k],x[j]) for k in range(3))-lap(q[i,j])-s.diff(s.trace(q),x[i],x[j]))/2)
    xis=s.Matrix([fun('xi'+str(i)) for i in range(3)])
    pg=s.Matrix(3,3,lambda i,j:s.diff(xis[i],x[j])+s.diff(xis[j],x[i]))
    check('linearized Ricci pure spatial gauge', Ric(pg))
    check('linearized stationary mass equation', Ric(-2*n*s.eye(3))-s.hessian(n,x)-s.eye(3)*lap(n))
    K=s.Matrix(3,3,lambda i,j:(s.diff(bv[i],x[j])+s.diff(bv[j],x[i]))/2)
    check('linearized momentum operator', s.Matrix([sum(s.diff(K[j,i],x[j]) for j in range(3))-s.diff(s.trace(K),x[i]) for i in range(3)])-(s.Matrix([lap(q) for q in bv])-grad(div(bv)))/2)
    chi=fun('chi')
    check('stationary time gauge shift is momentum pure gauge', s.Matrix([lap(q) for q in grad(chi)])-grad(lap(chi)))
    ch=fun('chi_phase')
    zp=s.exp(-s.I*om*ch)*z
    check('time gauge covariant profile derivative', s.diff(zp,x[0])-s.exp(-s.I*om*ch)*(s.diff(z,x[0])-s.I*om*s.diff(ch,x[0])*z))
    # Full contravariant metric under t'=t+chi(x), with stationary x unchanged.
    cv=s.Matrix(s.symbols('chi0:3',real=True)); jt=s.eye(4); jt[0,1:4]=cv.T
    gp=jt*gi*jt.T
    check('exact stationary time coordinate inverse metric', gp[0,0]-gi[0,0]-2*gi[0,1:4].dot(cv)-cv.dot(gi[1:4,1:4]*cv))
    phi=fun('phi'); uu=s.exp(s.I*phi)*u
    check('full PDE phase conjugation retains Coulomb term', s.exp(-s.I*phi)*(lap(uu)+om**2*uu)-lap(u)-2*s.I*grad(phi).dot(grad(u))-(om**2+s.I*lap(phi)-grad(phi).dot(grad(phi)))*u)
    # Exact spherical reduction, not a theorem for the general tensor density.
    nr, hr, fr, Br, yr = [s.Function(q)(R) for q in ('N','H','f','B','y')]
    zr=s.Function('z')(R); ell=s.symbols('ell', integer=True)
    Cr=nr*hr*fr
    raw=s.diff(Cr*R**2*s.diff(zr,R),R)/(nr*hr**3*fr*R**2)-ell*(ell+1)*zr/(hr**2*R**2)+om**2*zr/nr**2
    target=s.diff(zr,R,2)+(2/R+s.diff(Cr,R)/Cr)*s.diff(zr,R)+(om**2*hr**2/nr**2-ell*(ell+1)/R**2)*zr
    check('exact spherical Stage G operator', hr**2*raw-target)
    zz=s.exp(-Br/2)*yr/R
    raw2=s.diff(zz,R,2)+(2/R+s.diff(Br,R))*s.diff(zz,R)+(om**2*hr**2/nr**2-ell*(ell+1)/R**2)*zz
    target2=s.diff(yr,R,2)+(om**2*hr**2/nr**2-ell*(ell+1)/R**2-s.diff(Br,R,2)/2-s.diff(Br,R)/R-s.diff(Br,R)**2/4)*yr
    check('exact radial density conjugation', R*s.exp(Br/2)*raw2-target2)
    k,c=s.symbols('k c', real=True, nonzero=True)
    th=k*R+c*s.log(R)/(2*k)
    for sign in (1,-1):
        wave=s.exp(sign*s.I*th)
        check('Coulomb logarithmic phase residual '+str(sign), (s.diff(wave,R,2)+(k*k+c/R)*wave)/wave+c*c/(4*k*k*R*R)+sign*s.I*c/(2*k*R*R))
    # Pointwise original-stress Komar / Noether cancellation, arbitrary shift.
    rr=s.Matrix(s.symbols('r0:3',real=True)); dx=s.Matrix(s.symbols('u0:3',real=True)); dy=s.Matrix(s.symbols('v0:3',real=True))
    u0,v0,V=s.symbols('u v V',real=True); zz=u0+s.I*v0; dz=dx+s.I*dy
    p=(s.I*om*zz-b.dot(dz))/N; q=-b.dot(rr)/N
    U=q*q+f*s.expand(s.conjugate(p)*p)
    j=-q*rr-s.Matrix([f*s.re(s.conjugate(p)*dd) for dd in dz])
    Q=f*s.im(s.conjugate(zz)*p)
    check('general shift Komar charge cancellation', N*U-b.dot(j)-om*Q)
    X=s.symbols('X',real=True); eps=(U+X)/2+V; S=(3*U-X)/2-3*V
    check('trace reversed Komar integrand', N*(eps+S)-2*b.dot(j)-(2*om*Q-2*N*V))
    check('normal Noether density from contravariant time current', N*f*s.im(s.conjugate(zz)*(-s.I*om*zz+b.dot(dz))/N**2)+Q)
    Ksym=s.Matrix(3,3,lambda i,j:fun(f'K{min(i,j)}{max(i,j)}')); kj=s.trace(Ksym)
    # Flat chart product rule is also covariant with metric compatible D.
    A=Ksym-s.eye(3)*kj
    source=s.Matrix([sum(s.diff(A[j,i],x[j]) for j in range(3)) for i in range(3)])
    check('momentum contracted shift boundary product', div(A*bv)-bv.dot(source)-sum(A[j,i]*s.diff(bv[i],x[j]) for i in range(3) for j in range(3)))
    # Nontrivial version with arbitrary boundary normal and tensor.
    sv=s.Matrix(s.symbols('s0:3',real=True)); kt=s.Symbol('Ktrace')
    check('Komar surface from two ADM integrations', (b.T*(Ksym-s.eye(3)*kt)*sv)[0]+kt*b.dot(sv)-(b.T*Ksym*sv)[0])
    check('phase real energy boundary product', vb*div(aa*grad(u))-div(vb*aa*grad(u))+grad(vb).dot(aa*grad(u)))
    # One-dimensional local dilation commutator; no spectral estimate follows.
    xx=x[0]; av=fun('comm_A'); bv1=fun('comm_B'); vv=fun('comm_V')
    def hop(q):
        return -s.diff(av*s.diff(q,xx),xx)-s.I*om*(2*bv1*s.diff(q,xx)+s.diff(bv1,xx)*q)+vv*q
    def dil(q): return xx*s.diff(q,xx)+q/2
    avec=2*av-xx*s.diff(av,xx); bvec=bv1-xx*s.diff(bv1,xx)
    expected=-s.diff(avec*s.diff(u,xx),xx)-s.I*om*(2*bvec*s.diff(u,xx)+s.diff(bvec,xx)*u)-xx*s.diff(vv,xx)*u
    check('local dilation commutator including shift',hop(dil(u))-dil(hop(u))-expected)
    print(f'PASS: Stage H {COUNT} exact local checks; analytical and unresolved claims are not symbolic tests.',flush=True)

if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('FAIL: Stage H execution:',type(exc).__name__,str(exc),flush=True)
        sys.exit(1)
