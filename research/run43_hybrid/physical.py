"""Hybridized quadratic radial SBP elements with a single-valued skeleton displacement."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import eigsh, splu, spsolve

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'research'))
from run43_consistent_momentum import Grid, G, scalar_flux_laplacian
from validate_toroidal_helical_dee_source import R, RMIN, RMAX
from diagnose_toroidal_boundary_convergence import Audit


def periodic(n, offset):
    rows = np.arange(n)
    return sp.csr_matrix((np.ones(n), (rows, (rows + offset) % n)), shape=(n,n))


class Hybrid:
    def __init__(self, n, conforming=False):
        self.n = n
        self.h = h = (RMAX-RMIN)/n
        dt = 2*np.pi/n
        ne = n//2
        # Broken endpoint, midpoint, endpoint; adjacent elements own distinct values.
        levels = (2*np.arange(ne)[:,None] + np.arange(3)).ravel()
        self.level = level = (levels[:,None]*n + np.arange(n)).ravel()
        self.size = size = len(level)
        self.r = r = RMIN + h*(level//n)
        theta = dt*(level % n)
        s, c = np.sin(theta), np.cos(theta)
        q = R+r*c
        self.q = q
        self.s, self.c = s, c
        local_d = np.array([[-3,4,-1],[-1,0,1],[1,-4,3]])/(2*h)
        weights = h*np.array([1,4,1])/3
        hd = weights[:,None]*local_d
        self.sbp = float(abs(hd+hd.T-np.diag([-1,0,1])).max())
        dr = sp.kron(sp.eye(ne),sp.kron(local_d,sp.eye(n)),format='csr')
        dt_op = sp.kron(sp.eye(3*ne),(periodic(n,1)-periodic(n,-1))/(2*dt))
        self.mass = mass = np.tile(np.repeat(weights,n),ne)*dt*r*q
        diag = sp.diags
        # At phi=0, Cartesian component ordering is cylindrical (R,phi,z).
        # At each other azimuth both field and basis rotate by the same angle.
        gx = diag(c)@dr-diag(s/r)@dt_op
        gz = diag(s)@dr+diag(c/r)@dt_op
        spin = diag(np.sin(dt)/(dt*q))
        zero = sp.csr_matrix((size,size))
        blocks = [[gx,zero,zero],[zero,gx,zero],[zero,zero,gx],
                  [zero,-spin,zero],[spin,zero,zero],[zero,zero,zero],
                  [gz,zero,zero],[zero,gz,zero],[zero,zero,gz]]
        gradient = [sp.hstack(b,format='csr') for b in blocks]
        div = gradient[0]+gradient[4]+gradient[8]
        strain = [[gradient[3*j+i]+gradient[3*i+j] for j in range(3)] for i in range(3)]
        for i in range(3):
            strain[i][i] -= (2/3)*div
        self.strain = L = sp.vstack(sum(strain,[]),format='csr')
        stiffness = .5*(L.T@diag(np.tile(mass,9))@L)
        # Energy form of angular fourth-difference stabilization.
        d2 = sp.kron(sp.eye(3*ne),(periodic(n,1)+periodic(n,-1)-2*sp.eye(n))/dt)
        damp_theta = .25*(d2.T@diag(mass/r**2)@d2)
        damp_phi = diag(.25*mass*((2*np.cos(dt)-2)/dt)**2/q**2)
        stiffness += sp.block_diag([damp_theta+damp_phi,damp_theta+damp_phi,damp_theta])
        faces = np.flatnonzero(np.repeat(np.tile([True,False,True],ne),n))
        cols = np.arange(len(faces))
        select = sp.csr_matrix((np.ones(len(faces)),(cols,faces)),shape=(len(faces),size))
        self.faces = faces
        sign = np.tile(np.repeat([-1,1],n),ne)
        traction = sp.vstack([select@(diag(c)@strain[i][0]+diag(s)@strain[i][2]) for i in range(3)])
        traction = diag(np.tile(sign,3))@traction
        # All skeleton faces are genuine DOFs. Only prescribed skeleton data are eliminated.
        nf = (ne-1)*n
        self.nf = nf
        sf = level[faces]//(2*n)*n + level[faces]%n - n
        interior = (sf>=0)&(sf<nf)
        skeleton = sp.csr_matrix((np.ones(interior.sum()),(cols[interior],sf[interior])),shape=(len(faces),nf))
        jump = sp.hstack([sp.block_diag([select]*3),-sp.block_diag([skeleton]*3)])
        self.jump = jump.tocsr()
        t = sp.hstack([traction,sp.csr_matrix((3*len(faces),3*nf))])
        surface = diag(np.tile(dt*(r*q)[faces],3))
        cross = t.T@surface@jump
        penalty = 8/weights[0]
        self.K = (sp.block_diag([stiffness,sp.csr_matrix((3*nf,3*nf))])
                  -cross-cross.T+penalty*(jump.T@surface@jump)).tocsc()
        self.prolong = None
        if conforming:
            # Identify matching traces algebraically, including both boundary conditions.
            # Every volume and skeleton row maps to the same continuous nodal space.
            ids = np.arange((n+1)*n)
            traces = ids.reshape(n+1,n)[2:-1:2].ravel()
            def interior_map(indices):
                keep = (indices>=n)&(indices<n*n)
                return sp.csr_matrix((np.ones(keep.sum()),(np.flatnonzero(keep),indices[keep]-n)),
                                     shape=(len(indices),(n-1)*n))
            volume_map = sp.block_diag([interior_map(level)]*3)
            face_map = sp.block_diag([interior_map(traces)]*3)
            self.prolong = sp.vstack([volume_map,face_map],format='csr')
            self.K = (self.prolong.T@self.K@self.prolong).tocsc()
        self.lu = splu(self.K)
        self.theta = theta

    def solve(self, force):
        rhs = np.r_[-(self.mass*force).ravel(),np.zeros(3*self.nf)]
        if self.prolong is not None:
            rhs = self.prolong.T@rhs
        solution = self.lu.solve(rhs)
        x = solution if self.prolong is None else self.prolong@solution
        self.last_solution = x.copy()
        w = x[:3*self.size].reshape(3,-1)
        nodal = np.zeros((3,self.n+1,self.n))
        nodal[:,1::2] = w.reshape(3,-1,3,self.n)[:,:,1]
        nodal[:,2:-1:2] = x[3*self.size:].reshape(3,-1,self.n)
        a2 = np.sum((self.strain@w.ravel()).reshape(9,-1)**2,axis=0)
        norm2 = np.bincount(self.level,weights=a2)/np.bincount(self.level)
        return nodal,norm2,dict(linear_l2=float(np.linalg.norm(self.K@solution-rhs)),
               rhs_l2=float(np.linalg.norm(rhs)),max_displacement_jump=float(abs(self.jump@x).max()))

    def manufactured(self):
        r,s,c,q = self.r,self.s,self.c,self.q
        # W=f(r) e_z has div W=f' sin(theta).
        f = np.sin(np.pi*(r-RMIN)/(RMAX-RMIN))
        k = np.pi/(RMAX-RMIN)
        fp = k*np.cos(k*(r-RMIN))
        fpp = -k*k*f
        force = np.array([c*s*(fpp-fp/r)/3, np.zeros_like(r),
                         fpp+(1/r+c/q)*fp+(s*s*fpp+c*c*fp/r)/3])
        w,_,info = self.solve(force)
        exact = np.sin(np.pi*np.arange(self.n+1)/self.n)
        err = w.copy()
        err[2] -= exact[:,None]
        for label,region in [('whole',slice(None)),('first_interior',1),('last_interior',-2)]:
            info[label+'_rms'] = float(np.sqrt(np.mean(np.sum(err[:,region]**2,axis=0))))
        # Separate truncation on the analytic field from differentiated solution error.
        phi = np.arange(self.n)*2*np.pi/self.n
        wr,wp,wz = w[:,:,:,None]
        full = np.array(np.broadcast_arrays(wr*np.cos(phi)-wp*np.sin(phi),
                        wr*np.sin(phi)+wp*np.cos(phi),wz))
        exact_field = np.zeros_like(full)
        exact_field[2] = exact[:,None,None]
        audit = Audit(full.shape[1:])
        er,et = audit.er,audit.et
        fp = k*np.cos(k*(audit.r-RMIN))
        fpp = -k*k*np.sin(k*(audit.r-RMIN))
        analytic = (er*np.sin(audit.t)*fpp+et*np.cos(audit.t)*fp/audit.r)/3
        analytic[2] += fpp+(1/audit.r+np.cos(audit.t)/audit.q)*fp
        def div_stress(field):
            stress = audit.tensor(field,2)
            return np.array([sum(audit.grad(stress[i,j],2)[j] for j in range(3)) for i in range(3)])
        truncation = div_stress(exact_field)-analytic
        differentiated_error = div_stress(full-exact_field)
        nodes = np.arange(self.n+1)
        band = (nodes>0)&(audit.r[:,0,0]<=.10+1e-12)
        jac = np.broadcast_to(audit.r*audit.q,full.shape[1:])[band]
        def rms(vec):
            return float(np.sqrt(np.sum(jac*np.sum(vec[:,band]**2,axis=0))/jac.sum()))
        info.update(inner_frozen_exact_field_defect=rms(truncation),
                    inner_differentiated_solution_error=rms(differentiated_error),
                    inner_manufactured_total_defect=rms(truncation+differentiated_error))
        return info


def run(n,out,conforming=False):
    a = Hybrid(n,conforming)
    mms = a.manufactured()
    # Scalar face-flux stencil and original physical matter projections.
    g = Grid(n+1,n,4,build_vector=False)
    sl = np.arange(0,g.size,4)
    flat = (n+1)*n
    grid = SimpleNamespace(shape=(n+1,n,1),size=flat,r=g.r[sl],q=g.q[sl],t=g.t[sl])
    lap = scalar_flux_laplacian(grid).tocsr()
    ii = np.arange(n,flat-n)
    u = np.ones(flat)
    spatial = g.spatial[sl]
    source = g.S[:,sl]
    def vector(u):
        return a.solve(8*np.pi*G*(source*u**6)[:,a.level])
    history = []
    for outer in range(30):
        w,a2,info = vector(u)
        def residual():
            energy = .5*g.temporal + .5*spatial*u**-4 + g.U
            return lap@(u-1)+a2*u**-7/8+2*np.pi*G*energy*u**5,energy
        for inner in range(12):
            res,e = residual()
            if abs(res[ii]).max()<1e-10:
                break
            d = -7*a2*u**-8/8 + 2*np.pi*G*(-2*spatial+5*e*u**4)
            step = spsolve(lap[ii][:,ii]+sp.diags(d[ii]),-res[ii])
            alpha = 1.
            while (u[ii]+alpha*step).min()<=0:
                alpha *= .5
            u[ii] += alpha*step
        else:
            raise RuntimeError('scalar Newton failed')
        new_w,a2,updated = vector(u)
        change = float(abs(w-new_w).max())
        hmax = float(abs(residual()[0][ii]).max())
        history.append(dict(initial=info,updated=updated,change=change,Hmax=hmax,newton_steps=inner))
        if change<1e-12 and hmax<1e-9:
            break
    else:
        raise RuntimeError('coupled solver failed')
    phi = np.arange(n)*2*np.pi/n
    wr,wp,wz = new_w[:,:,:,None]
    field = np.array(np.broadcast_arrays(wr*np.cos(phi)-wp*np.sin(phi),
                      wr*np.sin(phi)+wp*np.cos(phi),wz))
    psi = np.broadcast_to(u.reshape(n+1,n,1),(n+1,n,n))
    path = out/f'fields_{n}.npz'
    np.savez_compressed(path,psi=psi,W_cartesian=field)
    np.savez_compressed(out/f'element_{n}.npz',broken_and_trace=a.last_solution,level=a.level)
    skew = a.K-a.K.T
    diagonal = sp.diags(1/np.sqrt(a.K.diagonal()))
    scaled = (diagonal@a.K@diagonal).tocsc()
    eig = eigsh(scaled,k=1,sigma=0,which='LM',return_eigenvectors=False,v0=np.ones(scaled.shape[0]))[0]
    row = dict(n=n,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),history=history,
               manufactured=mms,radial_sbp_defect=a.sbp,
               max_stiffness_skew=float(abs(skew).max()),scaled_eigenvalue_nearest_zero=float(eig),
               audit=Audit(psi.shape).summarize(psi,field))
    (out/f'result_{n}.json').write_text(json.dumps(row,indent=2)+'\n')
    print(n,mms,history[-1],flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--conforming',action='store_true',help='use identical elements with shared strong traces')
    args = parser.parse_args()
    args.directory.mkdir(parents=True,exist_ok=False)
    for n in (8,16,32,64):
        run(n,args.directory,args.conforming)
