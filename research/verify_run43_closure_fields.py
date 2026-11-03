"""Separate implementation of all twelve ORIGINAL frozen H/M2 checks.

Does not import Audit or solver derivative/tensor implementations. Reconstructs
original stress projections and evaluates explicit NumPy differences on arrays.
Checks against the separate frozen-code audit without changing its criterion.
"""
import sys,json,hashlib,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reproduce'))
from dee_stress_energy_3p1 import coeffs,potential
from validate_toroidal_helical_dee_source import R,RMIN,RMAX,M,RHO,PI,KTH,KPS


def evaluate(u,w):
    nr,n,_=u.shape;h=(RMAX-RMIN)/(nr-1);ha=2*np.pi/n
    r=np.linspace(RMIN,RMAX,nr)[:,None,None]
    t=np.arange(n)[None,:,None]*ha;p=np.arange(n)[None,None,:]*ha;q=R+r*np.cos(t)
    bc=lambda *x:np.array(np.broadcast_arrays(*x))
    er=bc(np.cos(t)*np.cos(p)+0*r,np.cos(t)*np.sin(p)+0*r,np.sin(t)+0*p+0*r)
    et=bc(-np.sin(t)*np.cos(p)+0*r,-np.sin(t)*np.sin(p)+0*r,np.cos(t)+0*p+0*r)
    ep=bc(-np.sin(p)+0*r+0*t,np.cos(p)+0*r+0*t,0*r+0*t+0*p)
    def grad(v):
        dr=np.empty_like(v);dr[1:-1]=(v[2:]-v[:-2])/(2*h)
        dr[0]=(-3*v[0]+4*v[1]-v[2])/(2*h)
        dr[-1]=(3*v[-1]-4*v[-2]+v[-3])/(2*h)
        return er*dr+et*(np.roll(v,-1,1)-np.roll(v,1,1))/(2*ha*r)+ep*(np.roll(v,-1,2)-np.roll(v,1,2))/(2*ha*q)
    d=np.array([grad(v) for v in w]);div=np.trace(d,axis1=0,axis2=1)
    a=d+d.swapaxes(0,1)-np.eye(3)[:,:,None,None,None]*(2/3)*div
    kt,kp,c=coeffs(RHO);kin=np.array([[1,0,0],[0,kt,c],[0,c,kp]])
    derivatives=np.zeros((3,3)+u.shape);direction=et/r-M*ep/q
    derivatives[1]=KTH*direction;derivatives[2]=KPS*direction
    energy=.5*float(PI@kin@PI)+potential(RHO)+.5*np.einsum('ab,ai...,bi...->...',kin,derivatives,derivatives)*u**-4
    source=-np.einsum('ab,a,bi...->i...',kin,PI,derivatives)
    momentum=u**-10*np.array([sum(grad(a[i,j])[j] for j in range(3)) for i in range(3)])-8*np.pi*.001*u**-4*source
    rf=(r[1:]+r[:-1])/2
    radial_flux=rf*(R+rf*np.cos(t))*(u[1:]-u[:-1])/h
    lap=(radial_flux[1:]-radial_flux[:-1])/(h*r[1:-1]*q[1:-1])
    angular_flux=(R+r*np.cos(t+ha/2))*(np.roll(u,-1,1)-u)/(ha*r)
    lap+=((angular_flux-np.roll(angular_flux,1,1))/(ha*r*q))[1:-1]
    lap+=((np.roll(u,-1,2)-2*u+np.roll(u,1,2))/(ha*ha*q*q))[1:-1]
    ham=-8*u[1:-1]**-5*lap-u[1:-1]**-12*np.sum(a*a,axis=(0,1))[1:-1]-16*np.pi*.001*energy[1:-1]
    weight=np.broadcast_to(r*q,u.shape)*u**6;rows={}
    for band,lo,hi in [('inner_boundary_band',.05,.10),('outer_boundary_band',.20,.25)]:
        mask=np.broadcast_to((r[1:-1]>=lo-1e-12)&(r[1:-1]<=hi+1e-12),u[1:-1].shape)
        wt=weight[1:-1][mask];den=np.sum(wt)
        rows[band]={'independent_flux_H':float(np.sqrt(np.sum(ham[mask]**2*wt)/den)),
                    'independent_tensor_M_order2':float(np.sqrt(np.sum(np.sum(momentum[:,1:-1]**2,axis=0)[mask]*wt)/den))}
    return rows


def main():
    directory=Path(sys.argv[1]);rows=[]
    frozen=json.loads((directory/'independent-audit.json').read_text())
    for k,n in enumerate((8,16,32,64)):
        path=directory/f'fields_{n}.npz'
        with np.load(path) as z:u,w=z['psi'],z['W_cartesian']
        result=evaluate(u,w)
        for band,metrics in result.items():
            for metric,value in metrics.items():
                expected=frozen['rows'][k]['audit'][band][metric]['volume_weighted_rms']
                np.testing.assert_allclose(value,expected,rtol=2e-7,atol=1e-12)
        rows.append(dict(n=n,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),norms=result))
    checks=[]
    for band in rows[0]['norms']:
        for metric in rows[0]['norms'][band]:
            for a,b in zip(rows,rows[1:]):
                p=math.log2(a['norms'][band][metric]/b['norms'][band][metric])
                checks.append(dict(band=band,metric=metric,pair=[a['n'],b['n']],order=p,passed=math.isfinite(p) and p>=1.5))
    report=dict(rows=rows,checks=checks,passed=sum(x['passed'] for x in checks),total=12,status='PASS' if all(x['passed'] for x in checks) else 'FAIL',criterion='All twelve finite orders >=1.5; unchanged frozen bands and volume weights',method='Independent explicit-array implementation, compared to frozen-code audit')
    (directory/'independent-implementation-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
    return 0 if report['status']=='PASS' else 1


if __name__=='__main__':sys.exit(main())
