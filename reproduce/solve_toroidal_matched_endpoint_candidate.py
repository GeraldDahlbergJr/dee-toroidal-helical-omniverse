"""Experimental endpoint closure only; frozen Run42 production solver is intact."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import splu
from solve_toroidal_helical_coupled_constraints import Grid,G
from investigate_toroidal_inner_operator import endpoint_cubic_matrix
from diagnose_toroidal_boundary_convergence import Audit
from validate_toroidal_helical_dee_source import RMIN,RMAX

class CandidateGrid(Grid):
    def __init__(self,nr,nt,nph):
        super().__init__(nr,nt,nph)
        dr=(RMAX-RMIN)/(nr-1)
        oldrad=self.d[0]
        newrad=sp.kron(sp.kron(sp.csr_matrix(endpoint_cubic_matrix(nr,dr,matched=True)),sp.eye(nt),format='csr'),sp.eye(nph),format='csr')
        # Only the first-derivative endpoint rows change. Direct second
        # derivatives, matter, geometry, boundary values and tolerances persist.
        self.d[0]=newrad
        self.lap=(self.lap+sp.diags(1/self.r+np.cos(self.t)/self.q)@(newrad-oldrad)).tocsr()
        er,et,ep=self.basis
        self.grad=[(sp.diags(er[i])@self.d[0]+sp.diags(et[i]/self.r)@self.d[1]+sp.diags(ep[i]/self.q)@self.d[2]).tocsr() for i in range(3)]
        L=self.lap[self.idx]@self.inject
        self.L=L.tocsc();self.lu=splu(self.L)
        self.vector=sp.bmat([[(L if i==j else sp.csr_matrix(L.shape))+(self.grad[i]@self.grad[j])[self.idx]@self.inject/3 for j in range(3)] for i in range(3)],format='csr')


def main():
    root=Path(__file__).parent;rows=[];hashes={}
    for n in (8,16,32):
        print(f'Candidate solve {n+1} x {n} x {n}',flush=True)
        grid=CandidateGrid(n+1,n,n);u,w,iterations=grid.solve()
        scalar=u.reshape(grid.shape);vector=w.reshape((3,)+grid.shape)
        a=Audit(grid.shape);regions=a.summarize(scalar,vector)
        # Supplementary matched endpoint audit cannot replace the original gate.
        a.drad[2]=endpoint_cubic_matrix(n+1,a.dr,matched=True)
        matched=a.summarize(scalar,vector)
        row={'resolution':list(grid.shape),'solver_diagnostics':grid.diagnostics(u,w),
             'coupled_iterations':iterations,'original_independent_audit':regions,
             'matched_endpoint_audit':matched}
        rows.append(row)
        name=f'toroidal_candidate_fields_{n+1}_{n}.npz';path=root/name
        np.savez_compressed(path,psi=scalar,W_cartesian=vector,resolution=grid.shape)
        hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        print(json.dumps({'resolution':row['resolution'],'expanded_max':row['solver_diagnostics']['expanded_momentum_solver_max_abs']}),flush=True)
    orders={}
    for route in ('original_independent_audit','matched_endpoint_audit'):
        orders[route]={}
        for region in ('inner_boundary_band','bulk','outer_boundary_band'):
            orders[route][region]={}
            for metric in ('independent_flux_H','independent_tensor_M_order2','independent_tensor_M_order4'):
                vals=[r[route][region][metric]['volume_weighted_rms'] for r in rows]
                orders[route][region][metric]=[float(np.log2(a/b)) for a,b in zip(vals,vals[1:])]
    gate={region:all(p>=1.5 for metric in ('independent_flux_H','independent_tensor_M_order2') for p in orders['original_independent_audit'][region][metric]) for region in ('inner_boundary_band','outer_boundary_band')}
    out={'status':'CANDIDATE_BOUNDARY_GATE_PASS_VALIDATION_PENDING' if all(gate.values()) else 'CANDIDATE_BOUNDARY_GATE_NOT_MET',
         'candidate':'matched leading truncation error at both radial first-derivative endpoints',
         'equations_and_seed':'same Run42 geometry, G, matter, boundary values, scalar/vector equations and tolerances',
         'production_solver_changed':False,'rows':rows,'observed_orders':orders,'original_boundary_gate':gate,
         'raw_field_sha256':hashes,'candidate_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'boundary_criterion':'unchanged: independent second-order H and M RMS order >=1.5 on both refinement pairs in both fixed boundary bands',
         'limitations':['Experimental discretization; no freeze, tag or production promotion.','Matched endpoint audit is supplementary and does not redefine the gate.','Axisymmetric stress on a finite shell remains; no generic 3D or evolution validation.']}
    (root/'toroidal_matched_endpoint_candidate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'gate':gate,'orders':orders},indent=2),flush=True)

if __name__=='__main__':main()
