"""Short 3-D periodic ADM Einstein--DEE evolution of smooth helical ring data.

General six-component spatial metric and extrinsic curvature; no imposed
symmetry during evolution. Unit lapse, zero shift, short-time gate only.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from toroidal_helical_initial_data import solve,deriv,L,G,A,B,C,LV,V

FINAL=.12; CFL=.20

def geometric(g,k,h):
    inv=np.linalg.inv(g)
    dg=np.stack([deriv(g,h,a) for a in range(3)],axis=-3) # ... a,i,j
    ch=np.empty_like(dg) # ... k,i,j
    for i in range(3):
        for j in range(3):
            tmp=dg[...,i,:,j]+dg[...,j,:,i]-dg[...,:,i,j]
            ch[..., :,i,j]=.5*np.einsum('...kl,...l->...k',inv,tmp)
    dch=np.stack([deriv(ch,h,a) for a in range(3)],axis=-4) # ... a,k,i,j
    trace=np.einsum('...kik->...i',ch)
    ric=np.einsum('...kkij->...ij',dch)
    for j in range(3): ric[..., :,j]-=deriv(trace,h,j)
    ric+=np.einsum('...kij,...k->...ij',ch,trace)
    ric-=np.einsum('...lik,...kjl->...ij',ch,ch)
    kup=np.einsum('...ij,...jk->...ik',inv,k)
    kt=np.trace(kup,axis1=-2,axis2=-1)
    scalar=np.einsum('...ij,...ij->...',inv,ric)
    ham=scalar+kt**2-np.einsum('...ij,...ji->...',kup,kup)
    mom=np.zeros(g.shape[:-2]+(3,))
    for j in range(3): mom+=deriv(kup[...,j,:],h,j)
    # Gamma^j_jl K^l_i - Gamma^l_ji K^j_l - partial_i K.
    mom+=np.einsum('...l,...li->...i',trace,kup)
    mom-=np.einsum('...lji,...jl->...i',ch,kup)
    for i in range(3): mom[...,i]-=deriv(kt,h,i)
    return inv,ric,ham,mom,kt

def matter(g,inv,f,p,h):
    grad=np.stack([deriv(f,h,a) for a in range(3)],axis=-1) # ... field,i
    rho=f[...,0]; kt=1+A*rho*rho; kp=1+B*rho*rho; mix=C*rho*rho
    fm=np.zeros(f.shape[:-1]+(3,3));fm[...,0,0]=1
    fm[...,1,1]=kt;fm[...,2,2]=kp;fm[...,1,2]=fm[...,2,1]=mix
    qpi=np.einsum('...ab,...a,...b->...',fm,p,p)
    qgrad=np.einsum('...ab,...ai,...ij,...bj->...',fm,grad,inv,grad)
    pot=LV/4*(rho*rho-V*V)**2
    e=.5*(qpi+qgrad)+pot
    lag=.5*(qpi-qgrad)-pot
    s=np.einsum('...ab,...a,...bi->...i',fm,p,grad)
    stress=np.einsum('...ab,...ai,...bj->...ij',fm,grad,grad)+lag[...,None,None]*g
    return e,s,stress,grad,fm

def rhs(state,h):
    g,k,f,p=state
    if np.min(np.linalg.eigvalsh(g))<=0: raise FloatingPointError('metric signature lost')
    inv,ric,_,_,kt=geometric(g,k,h)
    e,s,stress,grad,fm=matter(g,inv,f,p,h)
    tr=np.einsum('...ij,...ij->...',inv,stress)
    dk=ric+kt[...,None,None]*k-2*np.einsum('...ik,...kl,...lj->...ij',k,inv,k)
    dk-=8*np.pi*G*(stress-.5*g*(tr-e)[...,None,None])
    rho=f[...,0]; root=np.sqrt(np.linalg.det(g))
    def divergence(q):
        out=np.zeros(q.shape[:-1]);
        for i in range(3): out+=deriv(root*q[...,i],h,i)
        return out/root
    gradup=np.einsum('...ij,...aj->...ai',inv,grad)
    laprho=divergence(gradup[...,0,:])
    q=np.einsum('...i,...ij,...j->...',grad[...,1,:],inv,grad[...,1,:])
    r=np.einsum('...i,...ij,...j->...',grad[...,2,:],inv,grad[...,2,:])
    cross=np.einsum('...i,...ij,...j->...',grad[...,1,:],inv,grad[...,2,:])
    charge=A*(q-p[...,1]**2)+B*(r-p[...,2]**2)+2*C*(cross-p[...,1]*p[...,2])
    dp0=kt*p[...,0]-laprho+LV*rho*(rho*rho-V*V)+rho*charge
    trans=fm[...,1:,1:]
    cur=np.einsum('...ab,...b->...a',trans,p[...,1:])
    flux=np.einsum('...ab,...bi->...ai',trans,gradup[...,1:,:])
    divflux=np.stack([divergence(flux[...,i,:]) for i in range(2)],axis=-1)
    fdot=-2*rho[...,None,None]*p[...,0,None,None]*np.array([[A,C],[C,B]])
    trhs=kt[...,None]*cur-divflux-np.einsum('...ab,...b->...a',fdot,p[...,1:])
    dpt=np.linalg.solve(trans,trhs[...,None])[...,0]
    return -2*k,dk,-p,np.concatenate((dp0[...,None],dpt),axis=-1)

def add(state,step,dt):return tuple(x+dt*y for x,y in zip(state,step))
def rk4(state,h,dt):
    a=rhs(state,h);b=rhs(add(state,a,dt/2),h)
    c=rhs(add(state,b,dt/2),h);d=rhs(add(state,c,dt),h)
    return tuple(x+dt*(u+2*v+2*w+z)/6 for x,u,v,w,z in zip(state,a,b,c,d))

def initial(n):
    f,psi,mean_k,record=solve(n)
    g=np.zeros((n,n,n,3,3));k=np.zeros_like(g)
    for i in range(3):g[...,i,i]=psi**4;k[...,i,i]=mean_k*psi**4/3
    return (g,k,f,np.zeros_like(f)),record

def diagnostics(state,h):
    g,k,f,p=state
    inv,ric,ham,mom,kt=geometric(g,k,h)
    e,s,stress,_,_=matter(g,inv,f,p,h)
    ham-=16*np.pi*G*e; mom-=8*np.pi*G*s
    rho=f[...,0]
    return dict(hamiltonian_l2=float(np.sqrt(np.mean(ham*ham))),
                momentum_l2=float(np.sqrt(np.mean(mom*mom))),
                hamiltonian_max=float(np.max(np.abs(ham))),
                momentum_max=float(np.max(np.abs(mom))),
                rho_min=float(np.min(rho)),metric_eigenvalue_min=float(np.min(np.linalg.eigvalsh(g))),
                transport_det_min=float(np.min((1+A*rho*rho)*(1+B*rho*rho)-(C*rho*rho)**2)),
                curvature_scalar_max_abs=float(np.max(np.abs(np.einsum('...ij,...ij->...',inv,ric)))))

def run(n):
    h=L/n;state,conformal=initial(n); start=diagnostics(state,h)
    steps=int(np.ceil(FINAL/(CFL*h)));dt=FINAL/steps
    for _ in range(steps):state=rk4(state,h,dt)
    return state,dict(n=n,steps=steps,dt=dt,conformal_initial=conformal,
                      initial=start,final=diagnostics(state,h))

def main(resolutions=(12,24,48),output='coupled_dee_ring_3d_checkpoint.json'):
    rows=[];states=[]
    if len(resolutions)!=3 or resolutions[1]!=2*resolutions[0] or resolutions[2]!=2*resolutions[1]:
        raise ValueError('resolutions must be three successively doubled grid sizes')
    for n in resolutions:
        state,row=run(n);states.append(state);rows.append(row)
    errors=[]
    for coarse,fine in zip(states,states[1:]):
        errors.append(float(max(np.max(np.abs(u-v[::2,::2,::2])) for u,v in zip(coarse,fine))))
    out=dict(scope='full spatial-tensor 3-D periodic ADM evolution of smooth helical-scalar ring; short fixed-gauge test',
             parameters=dict(G=G,L=L,t_end=FINAL,cfl=CFL),runs=rows,
             adjacent_resolution_max_errors=errors,
             observed_order=float(np.log2(errors[0]/errors[1])),
             gates=dict(hamiltonian_refines=all(rows[i+1]['final']['hamiltonian_l2']<rows[i]['final']['hamiltonian_l2'] for i in range(2)),
                        momentum_refines=all(rows[i+1]['final']['momentum_l2']<rows[i]['final']['momentum_l2'] for i in range(2)),
                        metric_positive=all(row['final']['metric_eigenvalue_min']>0 for row in rows),
                        finite_branch=all(row['final']['rho_min']>0 for row in rows)))
    Path(__file__).with_name(output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resolutions',nargs=3,type=int,default=(12,24,48))
    parser.add_argument('--output',default='coupled_dee_ring_3d_checkpoint.json')
    args=parser.parse_args()
    main(args.resolutions,args.output)
