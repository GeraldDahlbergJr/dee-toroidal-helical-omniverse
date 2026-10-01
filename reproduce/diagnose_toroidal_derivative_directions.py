"""Hybrid derivative audits isolate radial/angular effects; never redefine gates."""
from pathlib import Path
import hashlib,io,json,zipfile
import numpy as np
from diagnose_toroidal_boundary_convergence import Audit

class HybridAudit(Audit):
    def __init__(self,shape,radial_order,angular_order):
        super().__init__(shape);self.radial_order=radial_order;self.angular_order=angular_order
    def grad(self,u,order):
        ur=np.einsum('ij,jkl->ikl',self.drad[self.radial_order],u)
        ut=self.angular_derivative(u,1,self.dt,self.angular_order)
        up=self.angular_derivative(u,2,self.dp,self.angular_order)
        return self.er*ur+self.et*ut/self.r+self.ep*up/self.q

def main():
    root=Path(__file__).parent;out={'status':'SUPPLEMENTARY_DIRECTION_AUDIT_ONLY','families':{}}
    with zipfile.ZipFile(root/'frozen/run42/artifacts.zip') as z:
        frozen=json.loads(z.read('toroidal_helical_coupled_constraints.json'))
        candidate=json.loads((root/'toroidal_direct_hessian_candidate.json').read_text())
        for family,checkpoint,prefix in [('frozen_run42',frozen,'toroidal_coupled_fields'),('direct_hessian_candidate',candidate,'toroidal_direct_hessian_fields')]:
            rows=[]
            for old in checkpoint['rows']:
                nr,nt,np_=old['resolution'];name=f'{prefix}_{nr}_{nt}.npz'
                raw=z.read(name) if family=='frozen_run42' else (root/name).read_bytes()
                assert hashlib.sha256(raw).hexdigest()==checkpoint['raw_field_sha256'][name]
                with np.load(io.BytesIO(raw)) as f:u,w=f['psi'],f['W_cartesian']
                row={'resolution':[nr,nt,np_],'field_sha256':hashlib.sha256(raw).hexdigest(),'inner_band_momentum_rms':{}}
                for radial,angular in [(2,2),(4,2),(2,4),(4,4)]:
                    a=HybridAudit(u.shape,radial,angular)
                    r=a.r[:,0,0];mask=np.broadcast_to(((r>.05+1e-12)&(r<=.1+1e-12))[:,None,None],u.shape)
                    row['inner_band_momentum_rms'][f'radial{radial}_angular{angular}']=a.norms(a.momentum(u,w,2),u,mask)['volume_weighted_rms']
                rows.append(row)
            orders={key:[float(np.log2(a['inner_band_momentum_rms'][key]/b['inner_band_momentum_rms'][key])) for a,b in zip(rows,rows[1:])] for key in rows[0]['inner_band_momentum_rms']}
            out['families'][family]={'rows':rows,'orders':orders}
    out['diagnostic_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out['scope']='Only differentiation of hash-verified saved fields changes. This is not directional grid refinement, does not modify the fields or independent gate, and cannot establish a solver remedy.'
    (root/'toroidal_derivative_directions.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({family:data['orders'] for family,data in out['families'].items()},indent=2))
if __name__=='__main__':main()
