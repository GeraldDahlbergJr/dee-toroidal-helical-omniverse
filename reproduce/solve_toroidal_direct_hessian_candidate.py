"""Experimental direct Cartesian Hessian vector operator on the same shell.

Continuum equation and source match Run42. The discrete grad(div W) is formed
from chart Hessians with toroidal connection terms and direct second derivatives,
instead of composing discrete Cartesian first derivatives. No production change.
"""
from pathlib import Path
import hashlib,json,os
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator, gmres
from solve_toroidal_helical_coupled_constraints import Grid,derivative
from diagnose_toroidal_boundary_convergence import Audit
from validate_toroidal_helical_dee_source import RMIN,RMAX

class DirectHessianGrid(Grid):
    def __init__(self,nr,nt,nph):
        super().__init__(nr,nt,nph)
        kron=lambda a,b,c:sp.kron(sp.kron(a,b,format='csr'),c,format='csr')
        ir,it,ip=[sp.eye(n,format='csr') for n in self.shape]
        drr=kron(derivative(nr,(RMAX-RMIN)/(nr-1),second=True),it,ip)
        dtt=kron(ir,derivative(nt,2*np.pi/nt,True,True),ip)
        dpp=kron(ir,it,derivative(nph,2*np.pi/nph,True,True))
        dr,dt,dp=self.d;r,t,q=self.r,self.t,self.q;diag=sp.diags
        # Hessian components in the orthonormal (er, et, ep) frame.
        components=[drr,
            diag(1/r**2)@dtt+diag(1/r)@dr,
            diag(1/q**2)@dpp+diag(np.cos(t)/q)@dr-diag(np.sin(t)/(r*q))@dt,
            diag(1/r)@(dr@dt)-diag(1/r**2)@dt,
            diag(1/q)@(dr@dp)-diag(np.cos(t)/q**2)@dp,
            diag(1/(r*q))@(dt@dp)+diag(np.sin(t)/q**2)@dp]
        er,et,ep=self.basis
        # Keep the six chart-Hessian operators, but do not assemble nine Cartesian
        # Hessian matrices plus a 3x3 block matrix at fine resolution.  The 65^3
        # hosted runner was being reclaimed during that memory-heavy construction.
        # This LinearOperator applies the identical discrete operator matrix-free;
        # equations, stencils, source, tolerances and boundary conditions are unchanged.
        self._hessian_components=components
        self._cartesian_hessian_coefficients=[
            [[er[i]*er[j],et[i]*et[j],ep[i]*ep[j],
              er[i]*et[j]+et[i]*er[j],er[i]*ep[j]+ep[i]*er[j],
              et[i]*ep[j]+ep[i]*et[j]] for j in range(3)] for i in range(3)]
        self._vector_lap=self.lap[self.idx]@self.inject
        nint=len(self.idx)
        def vector_matvec(x):
            xin=np.asarray(x).reshape(3,nint)
            full=[self.inject@xin[j] for j in range(3)]
            # Apply each chart Hessian component once per input component.
            hc=[[op@full[j] for op in self._hessian_components] for j in range(3)]
            out=[]
            for i in range(3):
                yi=self._vector_lap@xin[i]
                for j in range(3):
                    coeff=self._cartesian_hessian_coefficients[i][j]
                    hij=sum(coeff[k]*hc[j][k] for k in range(6))
                    yi=yi+hij[self.idx]/3
                out.append(yi)
            return np.concatenate(out)
        self.vector=LinearOperator((3*nint,3*nint),matvec=vector_matvec,dtype=float)


    def solve_vector(self,u,x0=None):
        rhs=(8*np.pi*G*u[self.idx]**6*self.S[:,self.idx]).ravel(); n=len(self.idx)
        pre=LinearOperator((3*n,3*n),matvec=lambda x:np.concatenate([self.lu.solve(y) for y in x.reshape(3,n)]))
        x,info=gmres(self.vector,rhs,M=pre,x0=None if x0 is None else x0[:,self.idx].ravel(),rtol=1e-10,atol=1e-13,restart=60,maxiter=100)
        if info: raise RuntimeError(f'vector GMRES failed: {info}')
        w=np.zeros((3,self.size)); w[:,self.idx]=x.reshape(3,n); return w

    def solve_restartable(self, checkpoint_path, max_new_iterations=None):
        checkpoint_path=Path(checkpoint_path); u=np.ones(self.size); w=None; start=0
        if checkpoint_path.exists():
            z=np.load(checkpoint_path); u=z['u']; w=z['w']; start=int(z['next_iteration'])
            print(f'Resuming coupled solve at iteration {start}',flush=True)
        stop=30 if max_new_iterations is None else min(30,start+max_new_iterations)
        for iteration in range(start,stop):
            w=self.solve_vector(u,w); a=self.longitudinal(w); a2=np.sum(a*a,axis=(0,1))
            for _ in range(12):
                e=self.energy(u); f=self.lap@(u-1)+a2*u**-7/8+2*np.pi*G*e*u**5; residual=f[self.idx]
                if np.max(abs(residual))<1e-10: break
                de=-2*self.spatial*u**-5
                jac=self.L+sp.diags((-7*a2*u**-8/8+2*np.pi*G*(de*u**5+5*e*u**4))[self.idx])
                from scipy.sparse.linalg import splu
                step=splu(jac.tocsc()).solve(-residual); alpha=1.
                while np.min(u[self.idx]+alpha*step)<=0: alpha/=2
                u[self.idx]+=alpha*step
            else: raise RuntimeError('Hamiltonian Newton failed')
            wnew=self.solve_vector(u,w); change=np.max(abs(wnew-w)); w=wnew
            anew=self.longitudinal(w); a2=np.sum(anew*anew,axis=(0,1))
            h=(self.lap@(u-1)+a2*u**-7/8+2*np.pi*G*self.energy(u)*u**5)[self.idx]
            np.savez_compressed(checkpoint_path,u=u,w=w,next_iteration=iteration+1,change=change,hmax=np.max(abs(h)))
            print(json.dumps({'checkpoint_iteration':iteration+1,'change':float(change),'Hmax':float(np.max(abs(h)))}),flush=True)
            if change<1e-12 and np.max(abs(h))<1e-9: return u,w,iteration+1,True
        if max_new_iterations is not None:
            return u,w,stop,False
        raise RuntimeError('coupled iteration failed')

def manufactured_check(levels=(8,16)):
    rows=[];c=np.array([.3,-.4,.5])[:,None]
    for n in levels:
        g=DirectHessianGrid(n+1,n,n);r=g.r
        v=(r-RMIN)*(RMAX-r);vp=RMIN+RMAX-2*r;e=np.exp(8*(r-RMIN))
        f=v*e;fp=e*(vp+8*v);fpp=e*(-2+16*vp+64*v);w=c*f
        er,et,ep=g.basis
        exact=c*(fpp+(1/r+np.cos(g.t)/g.q)*fp)+(
            er*fpp*np.sum(er*c,axis=0)+et*fp/r*np.sum(et*c,axis=0)+ep*fp*np.cos(g.t)/g.q*np.sum(ep*c,axis=0))/3
        actual=(g.vector@w[:,g.idx].ravel()).reshape(3,-1)
        error=np.full((3,g.size),np.nan);error[:,g.idx]=actual-exact[:,g.idx]
        audit=Audit(g.shape);first=np.zeros(g.shape,dtype=bool);first[1]=True
        inner=np.broadcast_to(((audit.r[:,0,0]>RMIN+1e-12)&(audit.r[:,0,0]<=.1+1e-12))[:,None,None],g.shape)
        err=error.reshape((3,)+g.shape)
        rows.append({'resolution':list(g.shape),'inner_first_row_rms':audit.norms(err,np.ones(g.shape),first)['volume_weighted_rms'],
                     'inner_band_rms':audit.norms(err,np.ones(g.shape),inner)['volume_weighted_rms']})
    orders={key:[float(np.log2(a[key]/b[key])) for a,b in zip(rows,rows[1:])] for key in ('inner_first_row_rms','inner_band_rms')}
    return {'rows':rows,'orders':orders,'gate':all(p>=1.5 for vals in orders.values() for p in vals)}


def main():
    root=Path(__file__).parent;rows=[];hashes={}
    progress_path=root/'toroidal_direct_hessian_candidate_progress.json'
    def preserve_progress(active_resolution=None, completed=False):
        payload={'status':'COMPLETE' if completed else 'IN_PROGRESS',
                 'active_resolution':active_resolution,
                 'completed_resolutions':[r['resolution'] for r in rows],
                 'rows':rows,'raw_field_sha256':hashes,
                 'candidate_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        progress_path.write_text(json.dumps(payload,indent=2)+'\\n')
    # Main solve first; analytic operator tests are independently callable.
    levels=tuple(int(x) for x in os.environ.get('DEE_CANDIDATE_LEVELS','8,16,32,64').split(',') if x.strip())
    preserve_progress()
    for n in levels:
        preserve_progress([n+1,n,n])
        print(f'Direct Hessian solve {n+1} x {n} x {n}',flush=True)
        g=DirectHessianGrid(n+1,n,n); checkpoint=root/f'toroidal_direct_hessian_restart_{n+1}_{n}.npz'
        max_new=int(os.environ.get('DEE_MAX_NEW_ITERATIONS','0')) or None
        u,w,iterations,converged=g.solve_restartable(checkpoint,max_new)
        if not converged:
            print(json.dumps({'status':'CHECKPOINT_STAGE_COMPLETE','resolution':list(g.shape),'next_iteration':iterations}),flush=True)
            preserve_progress([n+1,n,n])
            return
        scalar=u.reshape(g.shape);vector=w.reshape((3,)+g.shape)
        audit=Audit(g.shape);regions=audit.summarize(scalar,vector)
        row={'resolution':list(g.shape),'coupled_iterations':iterations,
             'solver_diagnostics':g.diagnostics(u,w),'independent_audit':regions}
        rows.append(row);name=f'toroidal_direct_hessian_fields_{n+1}_{n}.npz';p=root/name
        np.savez_compressed(p,psi=scalar,W_cartesian=vector,resolution=g.shape)
        hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
        preserve_progress()
        print(json.dumps({'resolution':list(g.shape),'expanded_residual':row['solver_diagnostics']['expanded_momentum_solver_max_abs']}),flush=True)
    # A single-resolution job intentionally emits fields/diagnostics only; convergence
    # is assessed after combining it with the retained lower-resolution evidence.
    if len(rows) < 2:
        out={'status':'SINGLE_RESOLUTION_COMPLETE','candidate':'direct Cartesian Hessian from toroidal chart second derivatives and connection terms',
             'production_solver_changed':False,'equations_and_seed':'Run42 continuum equations, G, geometry, fields, boundaries and validated solver tolerances retained',
             'rows':rows,'raw_field_sha256':hashes,'candidate_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'limitations':['Single-resolution execution; combine with retained lower-resolution evidence before convergence assessment.','No production promotion, freeze or tag.']}
        (root/'toroidal_direct_hessian_candidate.json').write_text(json.dumps(out,indent=2)+'\\n')
        preserve_progress(completed=True)
        print(json.dumps({'status':out['status'],'completed_resolutions':out['rows'][0]['resolution'] if rows else []},indent=2),flush=True)
        return
    orders={}
    for region in ('inner_boundary_band','bulk','outer_boundary_band'):
        orders[region]={}
        for metric in ('independent_flux_H','independent_tensor_M_order2','independent_tensor_M_order4'):
            vals=[r['independent_audit'][region][metric]['volume_weighted_rms'] for r in rows]
            orders[region][metric]=[float(np.log2(a/b)) for a,b in zip(vals,vals[1:])]
    gate={region:all(p>=1.5 for metric in ('independent_flux_H','independent_tensor_M_order2') for p in orders[region][metric]) for region in ('inner_boundary_band','outer_boundary_band')}
    out={'status':'CANDIDATE_BOUNDARY_GATE_PASS_VALIDATION_PENDING' if all(gate.values()) else 'CANDIDATE_BOUNDARY_GATE_NOT_MET',
         'candidate':'direct Cartesian Hessian from toroidal chart second derivatives and connection terms',
         'production_solver_changed':False,'equations_and_seed':'Run42 continuum equations, G, geometry, fields, boundaries and validated solver tolerances retained',
         'rows':rows,'observed_orders':orders,'original_boundary_gate':gate,'raw_field_sha256':hashes,
         'candidate_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'coarse_grid_audit_retained':True,'fine_grid_assessment':'Report all refinement pairs separately; a coarse-pair failure is preserved rather than erased.','limitations':['Experimental discretization; no production promotion, freeze or tag.','Original independent second-order audit and order>=1.5 gate retained; fourth-order audit is supplementary.','Finite shell with axisymmetric stress; no general 3D or evolution validation.']}
    (root/'toroidal_direct_hessian_candidate.json').write_text(json.dumps(out,indent=2)+'\\n')
    preserve_progress(completed=True)
    print(json.dumps({'status':out['status'],'gate':gate,'orders':orders},indent=2),flush=True)

if __name__=='__main__':main()
