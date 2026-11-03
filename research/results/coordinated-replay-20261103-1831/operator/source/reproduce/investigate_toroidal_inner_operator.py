"""Manufactured and frozen-field investigation; no production solver changes.

A three-point endpoint first derivative is second order, but composing it
with an adjacent centered derivative can create a first-order second derivative.
The four-point endpoint closures below are audit-only comparisons; one matches
the leading centered truncation term instead of maximizing endpoint order.
"""
from pathlib import Path
import hashlib,io,json,zipfile
import numpy as np
from diagnose_toroidal_boundary_convergence import Audit,radial_matrix
from validate_toroidal_helical_dee_source import R,RMIN,RMAX


def endpoint_cubic_matrix(n,h,matched=False):
    # matched=True sets sum(weights*offset**3)=1, matching the centered
    # first derivative error +h**2*f'''/6 on BOTH radial endpoints.
    matrix=radial_matrix(n,h,2)
    for index,js in [(0,np.arange(4)),(n-1,np.arange(n-4,n))]:
        offsets=js-index
        weights=np.linalg.solve(np.array([offsets.astype(float)**k for k in range(4)]),np.array([0.,1.,0.,1. if matched else 0.]))/h
        matrix[index]=0.;matrix[index,js]=weights
    return matrix


def expanded(a,w):
    r,t,q=a.r,a.t,a.q
    dw=np.array([a.grad(v,2) for v in w]);div=np.einsum('ii...->...',dw)
    graddiv=a.grad(div,2)
    result=np.full(w.shape,np.nan)
    for i,v in enumerate(w):
        vr=np.einsum('ij,jkl->ikl',a.drad[2],v)
        lap=(v[2:]-2*v[1:-1]+v[:-2])/a.dr**2
        lap+=(1/r+np.cos(t)/q)[1:-1]*vr[1:-1]
        lap+=((np.roll(v,-1,1)-2*v+np.roll(v,1,1))/a.dt**2/r**2)[1:-1]
        lap+=(-np.sin(t)/(r*q)*a.angular_derivative(v,1,a.dt,2))[1:-1]
        lap+=((np.roll(v,-1,2)-2*v+np.roll(v,1,2))/a.dp**2/q**2)[1:-1]
        result[i,1:-1]=lap+graddiv[i,1:-1]/3
    return result


def tensor_div(a,w):
    tensor=a.tensor(w,2)
    return np.array([sum(a.grad(tensor[i,j],2)[j] for j in range(3)) for i in range(3)])


def region_norms(a,v,u):
    rv=a.r[:,0,0];out={}
    for name,radial in {
        'inner_first_row':np.arange(len(rv))==1,
        'inner_band':(rv>RMIN+1e-12)&(rv<=.1+1e-12),
        'bulk':(rv>=.1-1e-12)&(rv<=.2+1e-12),
        'outer_band':(rv>=.2-1e-12)&(rv<RMAX-1e-12),
    }.items():
        mask=np.broadcast_to(radial[:,None,None],a.shape)
        out[name]=a.norms(v,u,mask)['volume_weighted_rms']
    return out


def main():
    manufactured=[];radial=[];c=np.array([.3,-.4,.5])[:,None,None,None]
    for n in (8,16,32,64):
        a=Audit((n+1,n,n));r=a.r
        v=(r-RMIN)*(RMAX-r);vp=RMIN+RMAX-2*r
        exponential=np.exp(8*(r-RMIN));f=v*exponential
        fp=exponential*(vp+8*v);fpp=exponential*(-2+16*vp+64*v)
        w=np.broadcast_to(c*f,(3,)+a.shape).copy();u=np.ones(a.shape)
        erdot=np.sum(a.er*c,axis=0);etdot=np.sum(a.et*c,axis=0);epdot=np.sum(a.ep*c,axis=0)
        truth=c*(fpp+(1/r+np.cos(a.t)/a.q)*fp)+(
            a.er*fpp*erdot+a.et*fp/r*etdot+a.ep*fp*np.cos(a.t)/a.q*epdot)/3
        row={'resolution':list(a.shape),'standard_tensor_error':region_norms(a,tensor_div(a,w)-truth,u),
             'expanded_error':region_norms(a,expanded(a,w)-truth,u)}
        a.drad[2]=endpoint_cubic_matrix(n+1,a.dr)
        row['cubic_endpoint_tensor_error']=region_norms(a,tensor_div(a,w)-truth,u)
        a.drad[2]=endpoint_cubic_matrix(n+1,a.dr,matched=True)
        row['matched_endpoint_tensor_error']=region_norms(a,tensor_div(a,w)-truth,u)
        manufactured.append(row)
    # Isolate radial closure without angular errors; extend to asymptotic grids.
    for n in (9,17,33,65,129,257):
        x=np.linspace(RMIN,RMAX,n);h=x[1]-x[0];v=(x-RMIN)*(RMAX-x);vp=RMIN+RMAX-2*x
        f=v*np.exp(8*(x-RMIN));truth=np.exp(8*(x-RMIN))*(-2+16*vp+64*v)
        row={'nr':n}
        for name,d in [('standard',radial_matrix(n,h,2)),('cubic_endpoint',endpoint_cubic_matrix(n,h)),('matched_endpoint',endpoint_cubic_matrix(n,h,matched=True))]:
            error=d@d@f-truth
            row[name+'_first_inner_abs']=float(abs(error[1]))
        radial.append(row)
    root=Path(__file__).parent;archive=root/'frozen/run42/artifacts.zip'
    manifest=json.loads((root/'frozen/run42/manifest.json').read_text())
    assert 'sha256:'+hashlib.sha256(archive.read_bytes()).hexdigest()==manifest['source']['artifact_digest']
    frozen=[]
    with zipfile.ZipFile(archive) as z:
        checkpoint=json.loads(z.read('toroidal_helical_coupled_constraints.json'))
        for old in checkpoint['rows']:
            nr,nt,np_=old['resolution'];name=f'toroidal_coupled_fields_{nr}_{nt}.npz';raw=z.read(name)
            assert hashlib.sha256(raw).hexdigest()==checkpoint['raw_field_sha256'][name]
            with np.load(io.BytesIO(raw)) as fields:u,w=fields['psi'],fields['W_cartesian']
            a=Audit(u.shape);tensor=tensor_div(a,w);exp=expanded(a,w)
            dw=np.array([a.grad(v,2) for v in w]);div=np.einsum('ii...->...',dw)
            lap_direct=exp-a.grad(div,2)/3
            lap_composed=np.array([sum(a.grad(a.grad(v,2)[j],2)[j] for j in range(3)) for v in w])
            commutator=np.array([sum(a.grad(a.grad(w[j],2)[i],2)[j]-a.grad(a.grad(w[j],2)[j],2)[i] for j in range(3)) for i in range(3)])
            lap_defect=lap_composed-lap_direct
            np.testing.assert_allclose((tensor-exp)[:,1:-1],(lap_defect+commutator)[:,1:-1],atol=1e-11,rtol=1e-9)
            from solve_toroidal_helical_coupled_constraints import G
            physical=u**-10*tensor-8*np.pi*G*u**-4*a.S
            row={'resolution':list(u.shape),'field_sha256':hashlib.sha256(raw).hexdigest(),
                'standard_physical_M':region_norms(a,physical,u),
                'expanded_vs_tensor_defect':region_norms(a,u**-10*(tensor-exp),u),
                'composed_vs_direct_laplacian_defect':region_norms(a,u**-10*lap_defect,u),
                'cartesian_derivative_commutator':region_norms(a,u**-10*commutator,u)}
            a.drad[2]=endpoint_cubic_matrix(nr,a.dr)
            changed=u**-10*tensor_div(a,w)-8*np.pi*G*u**-4*a.S
            row['audit_only_cubic_endpoint_M']=region_norms(a,changed,u)
            a.drad[2]=endpoint_cubic_matrix(nr,a.dr,matched=True)
            changed=u**-10*tensor_div(a,w)-8*np.pi*G*u**-4*a.S
            row['audit_only_matched_endpoint_M']=region_norms(a,changed,u)
            frozen.append(row)
    def orders(rows,key,region):
        vals=[row[key][region] for row in rows]
        return [float(np.log2(a/b)) for a,b in zip(vals,vals[1:])]
    out={'status':'INVESTIGATION_ONLY_NO_SOLVER_CHANGE','source_run':42,
         'archive_digest':manifest['source']['artifact_digest'],
         'investigation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'manufactured_field':'W=c*(r-rmin)*(rmax-r)*exp(8*(r-rmin)); c=(0.3,-0.4,0.5); exact Cartesian vector Laplacian derived from Hessian(r)',
         'manufactured_rows':manufactured,'radial_closure_rows':radial,'frozen_rows':frozen,
         'manufactured_inner_orders':{key:orders(manufactured,key,'inner_first_row') for key in ['standard_tensor_error','expanded_error','cubic_endpoint_tensor_error','matched_endpoint_tensor_error']},
         'frozen_inner_band_orders':{key:orders(frozen,key,'inner_band') for key in ['standard_physical_M','expanded_vs_tensor_defect','audit_only_cubic_endpoint_M','audit_only_matched_endpoint_M']},
         'radial_first_inner_orders':{key:[float(np.log2(a[key]/b[key])) for a,b in zip(radial,radial[1:])] for key in ['standard_first_inner_abs','cubic_endpoint_first_inner_abs','matched_endpoint_first_inner_abs']},
         'scope':'Cubic endpoint comparison changes only the audit operator. Frozen solved fields and production equations, source, boundaries and solver remain unchanged. A modified solver and fresh convergence study would be required for a remedy.'}
    (root/'toroidal_inner_operator_investigation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','manufactured_inner_orders','frozen_inner_band_orders','radial_first_inner_orders']},indent=2))

if __name__=='__main__':main()
