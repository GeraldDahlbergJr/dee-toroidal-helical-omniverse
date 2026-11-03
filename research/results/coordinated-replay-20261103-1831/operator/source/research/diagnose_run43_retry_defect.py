"""Separate frozen endpoint truncation from independently measured residuals.

The matched audit below is supplementary only. The official twelve checks
continue to use the byte-identical original Audit with its original closure.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reproduce'))
from diagnose_toroidal_boundary_convergence import Audit
from investigate_toroidal_inner_operator import endpoint_cubic_matrix
from validate_original_physics import Independent


def main():
    root=Path(__file__).resolve().parents[1];rows=[]
    for variant,folder in [('preserved-12of12-baseline','continuation-20261103-fresh'),('radial-principal-retry','radial-retry-20261103-principal'),('fourth-vector-retry','radial-retry-20261103-fourth-vector')]:
        for n in (8,16,32,64):
            path=root/'research/results'/folder/f'fields_{n}.npz'
            with np.load(path) as z:u,w=z['psi'],z['W_cartesian']
            a=Audit(u.shape);frozen=a.momentum(u,w,2)
            b=Audit(u.shape);b.drad[2]=endpoint_cubic_matrix(n+1,b.dr,True)
            matched=b.momentum(u,w,2);endpoint_defect=frozen-matched
            ref=Independent(n);dw=np.array([ref.grad(v) for v in w]);div=np.trace(dw,axis1=0,axis2=1)
            tensor=dw+dw.swapaxes(0,1)-2/3*np.eye(3)[:,:,None,None,None]*div
            da=np.array([sum(ref.grad(tensor[i,j])[j] for j in range(3)) for i in range(3)])
            independent=u**-10*da-8*np.pi*.001*u**-4*a.S
            bands={}
            for label,lo,hi in [('inner',.05,.10),('outer',.20,.25)]:
                mask=np.broadcast_to((a.r>=lo-1e-12)&(a.r<=hi+1e-12),u.shape).copy();mask[[0,-1]]=False
                norm=lambda v:a.norms(v,u,mask)['volume_weighted_rms']
                weights=np.broadcast_to(a.r*a.q,u.shape)*u**6
                squares=np.sum(frozen*frozen,axis=0)*weights
                edge=np.zeros(u.shape,bool);edge[1 if label=='inner' else -2]=True
                bands[label]=dict(frozen_rms=norm(frozen),supplementary_matched_rms=norm(matched),endpoint_closure_difference_rms=norm(endpoint_defect),independent_high_order_rms=norm(independent),frozen_minus_high_order_rms=norm(frozen-independent),first_interior_row_fraction_of_frozen_squared_norm=float(np.sum(squares[edge&mask])/np.sum(squares[mask])))
            rows.append(dict(variant=variant,n=n,field_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bands=bands))
    report=dict(scope='Residual decomposition only. No substituted gate or acceptance threshold.',rows=rows)
    p=root/'research/results/radial-retry-20261103/residual-decomposition.json';p.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
