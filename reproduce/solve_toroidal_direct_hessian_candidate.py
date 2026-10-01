"""Experimental direct Cartesian Hessian vector operator on the same shell.

Continuum equation and source match Run42. The discrete grad(div W) is formed
from chart Hessians with toroidal connection terms and direct second derivatives,
instead of composing discrete Cartesian first derivatives. No production change.
"""
from pathlib import Path
import hashlib,json,os
import numpy as np
from scipy import sparse as sp
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
        self.hessian=[]
        for i in range(3):
            row=[]
            for j in range(3):
                coefficients=[er[i]*er[j],et[i]*et[j],ep[i]*ep[j],
                    er[i]*et[j]+et[i]*er[j],er[i]*ep[j]+ep[i]*er[j],et[i]*ep[j]+ep[i]*et[j]]
                row.append(sum(diag(c)@op for c,op in zip(coefficients,components)).tocsr())
            self.hessian.append(row)
        L=self.lap[self.idx]@self.inject
        self.vector=sp.bmat([[(L if i==j else sp.csr_matrix(L.shape))+self.hessian[i][j][self.idx]@self.inject/3 for j in range(3)] for i in range(3)],format='csr')


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
        g=DirectHessianGrid(n+1,n,n);u,w,iterations=g.solve()
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
