"""Continuous quartic radial Galerkin reference for the Run 43 equivariant sector."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import splu, spsolve

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "reproduce"))
from solve_toroidal_helical_coupled_constraints import G
from validate_toroidal_helical_dee_source import R, RMIN, RMAX, RHO, M, KTH, KPS, PI
from dee_stress_energy_3p1 import coeffs, potential


def quadrature(ne):
    x, weights = np.polynomial.legendre.leggauss(6)
    nodes = np.linspace(-1, 1, 5)
    width = (RMAX-RMIN)/ne
    values, derivatives = [], []
    for j, node in enumerate(nodes):
        others = np.delete(nodes,j)
        poly = np.polynomial.Polynomial.fromroots(others)
        poly = poly/poly(node)
        values.append(poly(x))
        derivatives.append(poly.deriv()(x)*2/width)
    columns = np.broadcast_to(4*np.arange(ne)[:,None,None]+np.arange(5),(ne,6,5))
    rows = np.broadcast_to(np.arange(6*ne).reshape(ne,6,1),columns.shape)
    def assemble(block):
        data = np.broadcast_to(np.array(block).T, columns.shape)
        return sp.csr_matrix((data.ravel(),(rows.ravel(),columns.ravel())),shape=(6*ne,4*ne+1))
    r = RMIN + width*(np.arange(ne)[:,None]+(x+1)/2)
    return assemble(values),assemble(derivatives),r.ravel(),np.tile(weights*width/2,ne)


class Reference:
    def __init__(self,n,nr):
        self.n, self.nr = n,nr
        h = 2*np.pi/n
        E,D,r,wr = quadrature(nr//4)
        eye = sp.eye(n)
        theta = np.arange(n)*h
        s,c = np.tile(np.sin(theta),len(r)),np.tile(np.cos(theta),len(r))
        self.r = r = np.repeat(r,n)
        self.c,self.s = c,s
        self.q = q = R+r*c
        self.mass = mass = np.repeat(wr,n)*h*r*q
        self.size = (nr+1)*n
        self.free = np.arange(n,nr*n)
        self.vfree = np.concatenate([self.free+i*self.size for i in range(3)])
        right = sp.csr_matrix((np.ones(n),(np.arange(n),(np.arange(n)+1)%n)),shape=(n,n))
        diff = (right-right.T)/(2*h)
        second = (right+right.T-2*eye)/h
        diag = sp.diags
        self.P = P = sp.kron(E,eye,format='csr')
        dr = sp.kron(D,eye,format='csr')
        dt = sp.kron(E,diff,format='csr')
        damp = sp.kron(E,second,format='csr')
        gx = diag(c)@dr-diag(s/r)@dt
        gz = diag(s)@dr+diag(c/r)@dt
        spin = diag(np.sin(h)/(q*h))@P
        zero = sp.csr_matrix(gx.shape)
        deriv = [sp.hstack(b,format='csr') for b in
                 [(gx,zero,zero),(zero,gx,zero),(zero,zero,gx),
                  (zero,-spin,zero),(spin,zero,zero),(zero,zero,zero),
                  (gz,zero,zero),(zero,gz,zero),(zero,zero,gz)]]
        div = deriv[0]+deriv[4]+deriv[8]
        components = []
        for i in range(3):
            for j in range(3):
                a = deriv[3*i+j]+deriv[3*j+i]
                if i==j:
                    a -= (2/3)*div
                components.append(a[:,self.vfree])
        self.L = sp.vstack(components,format='csr')
        angular = damp.T@diag(.25*mass/r**2)@damp
        azimuth = P.T@diag(.25*mass*((2*np.cos(h)-2)/h)**2/q**2)@P
        regularizer = sp.block_diag([angular+azimuth,angular+azimuth,angular],format='csr')
        self.K = (.5*self.L.T@diag(np.tile(mass,9))@self.L+regularizer[self.vfree][:,self.vfree]).tocsc()
        self.factor = splu(self.K)
        # Scalar theta face flux uses q at the angular midpoint. Its damping
        # multiplies the fourth difference outside the stencil, as in the hybrid.
        face = sp.kron(E,(right-eye)/h,format='csr')
        qface = R+r*(c*np.cos(h/2)-s*np.sin(h/2))
        mass_r2 = mass/r**2
        fourth = sp.kron(E,second@second,format='csr')
        scalar = dr.T@diag(mass)@dr + face.T@diag(mass_r2*qface/q)@face + P.T@diag(.25*mass_r2)@fourth
        self.A = scalar[self.free][:,self.free].tocsc()
        self.B = P[:,self.free]
        kt,kp,cross = coeffs(RHO)
        kin = np.array([[1,0,0],[0,kt,cross],[0,cross,kp]])
        gradients = np.zeros((3,3,len(r)))
        direction = np.array([-s/r,-M/q,c/r])
        gradients[1],gradients[2] = KTH*direction,KPS*direction
        self.temporal = PI@kin@PI
        self.spatial = np.einsum('ab,ain,bin->n',kin,gradients,gradients)
        self.source = -np.einsum('ab,a,bin->in',kin,PI,gradients)

    def vector(self,force):
        rhs = -np.asarray([self.B.T@(self.mass*v) for v in force]).ravel()
        x = self.factor.solve(rhs)
        w = np.zeros(3*self.size)
        w[self.vfree] = x
        norm2 = ((self.L@x).reshape(9,-1)**2).sum(axis=0)
        return w.reshape(3,self.nr+1,self.n),norm2,float(np.linalg.norm(self.K@x-rhs))

    def manufactured(self):
        # W=f e_z. The analytic radial flux derivative and exact angular
        # semidiscrete adjoint force keep angular truncation out of this test.
        r,s,c,q = self.r,self.s,self.c,self.q
        width = RMAX-RMIN
        x = (r-RMIN)/width
        f = lambda x: np.exp(x)*x*(1-x)
        fp = np.exp(x)*(1-x-x*x)/width
        fpp = np.exp(x)*(-3*x-x*x)/width**2
        # Add a second mode so f''(RMIN) is nonzero.
        fp += np.pi*np.cos(2*np.pi*x)/width + (1-2*x)/width
        fpp += -2*np.pi**2*np.sin(2*np.pi*x)/width**2 - 2/width**2
        er = np.array([c,np.zeros_like(c),s])
        et = np.array([-s,np.zeros_like(c),c])
        ez = np.array([[0],[0],[1]])
        radial = (fpp+(1/r+c/q)*fp)*(ez+er*s/3)
        flux = (q*fp*(er*c-et*(2*s/3))).reshape(3,-1,self.n)
        h = 2*np.pi/self.n
        angular = ((np.roll(flux,-1,axis=2)-np.roll(flux,1,axis=2))/(2*h)).reshape(3,-1)/(r*q)
        force = radial+angular
        force[0] += np.sin(h)/h*(2*s/3)*fp/q
        w,_,residual = self.vector(force)
        t = np.linspace(0,1,self.nr+1)
        exact = f(t)+.5*np.sin(2*np.pi*t)+t*(1-t)
        w[2] -= exact[:,None]
        at_audit = w[:,::self.nr//self.n]
        return dict(wall_curvature=-2/width**2,linear_residual=residual,
                    max_node_error=float(abs(w).max()),first_audit_row_error=float(abs(at_audit[:,1]).max()))

    def coupled(self):
        y = np.zeros(len(self.free))
        history = []
        def vector(y):
            return self.vector(8*np.pi*G*self.source*(1+self.B@y)**6)
        for outer in range(20):
            w,a2,linear = vector(y)
            def scalar():
                u = 1+self.B@y
                energy = .5*(self.temporal+self.spatial*u**-4)+potential(RHO)
                body = a2*u**-7/8+2*np.pi*G*energy*u**5
                residual = self.A@y-self.B.T@(self.mass*body)
                return u,energy,residual
            for it in range(10):
                u,e,res = scalar()
                slope = -7*a2*u**-8/8+2*np.pi*G*(5*e*u**4-2*self.spatial)
                step = spsolve(self.A-self.B.T@sp.diags(self.mass*slope)@self.B,-res)
                y += step
                if abs(step).max()<5e-15:
                    break
            else:
                raise RuntimeError('scalar Newton failed')
            new,a2,linear = vector(y)
            delta = float(abs(w-new).max())
            history.append(dict(iteration=outer+1,newton_steps=it+1,vector_l2=linear,
                                scalar_linf=float(abs(scalar()[2]).max()),change=delta))
            if delta<3e-16:
                nodal = np.ones(self.size)
                nodal[self.free] += y
                return nodal.reshape(self.nr+1,self.n),new,history
        raise RuntimeError('coupled iteration failed')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('--factors',type=int,nargs='+',default=[1,2,4])
    args=p.parse_args()
    args.directory.mkdir(exist_ok=True,parents=True)
    for n in (8,16,32,64):
        for factor in args.factors:
            nr=n*factor
            dest=args.directory/f'n{n}-r{nr}'
            dest.mkdir()
            ref=Reference(n,nr)
            mms=ref.manufactured()
            u,w,hist=ref.coupled()
            path=dest/'reduced.npz'
            np.savez_compressed(path,psi=u,W=w)
            record=dict(n=n,nr=nr,mms=mms,history=hist,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            (dest/'result.json').write_text(json.dumps(record,indent=2)+'\n')
            print(n,nr,mms,hist[-1],flush=True)
            del ref

if __name__=='__main__':
    main()
