"""Short plane-symmetric coupled Einstein--DEE ADM evolution (unit lapse, zero shift).

All six diagonal gamma_ij/K_ij components and three scalar fields evolve. The
1-D symmetry is a restriction of 3+1 GR, not a toroidal BSSN implementation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from nonlinear_initial_data import matter_projections

G=.005; L=8.; FINAL=.12; CFL=.12
A=.30; B=.20; C=.10; LV=4.; V=1.

def dx(f,h): return (np.roll(f,-1,axis=0)-np.roll(f,1,axis=0))/(2*h)
def dxx(f,h): return (np.roll(f,-1,axis=0)-2*f+np.roll(f,1,axis=0))/h**2

def geometry(g,k,h):
    ax=dx(g,h)/(2*g)
    gx,gy,gz=ax.T
    # Gamma^x_yy and Gamma^x_zz, plus Gamma^y_xy and Gamma^z_xz.
    cy=-dx(g[:,1],h)/(2*g[:,0]); cz=-dx(g[:,2],h)/(2*g[:,0])
    ry=dx(cy,h)+cy*(gx-gy+gz)
    rz=dx(cz,h)+cz*(gx+gy-gz)
    rx=-dx(gy+gz,h)+gx*(gy+gz)-gy**2-gz**2
    ric=np.column_stack((rx,ry,rz))
    scalar=np.sum(ric/g,axis=1)
    km=k/g; kt=np.sum(km,axis=1)
    ham=scalar+kt**2-np.sum(km**2,axis=1)
    mom=-dx(km[:,1]+km[:,2],h)+gy*(km[:,0]-km[:,1])+gz*(km[:,0]-km[:,2])
    return ric,ham,mom,kt

def matter(g,fields,pis,h):
    rho=fields[:,0]; grads=np.zeros((len(g),3,3)); grads[:,:,0]=dx(fields,h)
    metric=np.zeros((len(g),3,3)); inv=np.zeros_like(metric)
    for i in range(3): metric[:,i,i]=g[:,i]; inv[:,i,i]=1/g[:,i]
    return matter_projections(rho,pis[:,0],grads[:,0],pis[:,1],grads[:,1],
                              pis[:,2],grads[:,2],metric,inv,
                              lambda_v=LV,v=V,a=A,b=B,c=C)

def rhs(state,h):
    g,k,f,p=state
    if np.min(g)<=0: raise FloatingPointError('spatial metric lost positive signature')
    n=len(g); ric,_,_,kt=geometry(g,k,h)
    e,s,stress,trace=matter(g,f,p,h)
    sd=np.diagonal(stress,axis1=1,axis2=2)
    dg=-2*k
    dk=ric+kt[:,None]*k-2*k*k/g-8*np.pi*G*(sd-.5*g*(trace-e)[:,None])
    rho=f[:,0]; ft=1+A*rho**2; fp=1+B*rho**2; mix=C*rho**2
    kinetic=np.empty((n,2,2)); kinetic[:,0,0]=ft; kinetic[:,1,1]=fp
    kinetic[:,0,1]=kinetic[:,1,0]=mix
    if np.min(np.linalg.det(kinetic))<=0: raise FloatingPointError('transport kinetic matrix singular')
    root=np.sqrt(np.prod(g,axis=1)); flux=root/g[:,0]
    grad=dx(f,h)
    def div(q):
        if q.ndim == 1: return dx(flux*q,h)/root
        return dx(flux[:,None]*q,h)/root[:,None]
    drho=div(grad[:,0])
    # box rho - V' - (1/2) F'_AB (grad phi_A.grad phi_B)=0
    grad2=grad[:,1:] / np.sqrt(g[:,0,None])
    charge=A*(grad2[:,0]**2-p[:,1]**2)+B*(grad2[:,1]**2-p[:,2]**2)
    charge+=2*C*(grad2[:,0]*grad2[:,1]-p[:,1]*p[:,2])
    dp0=kt*p[:,0]-drho+LV*rho*(rho**2-V**2)+rho*charge
    current=np.einsum('nij,nj->ni',kinetic,p[:,1:])
    spatial=np.einsum('nij,nj->ni',kinetic,grad[:,1:])
    # evolve currents P_A=F_AB Pi_B; Fdot=-2 rho Pi_rho [[a,c],[c,b]]
    fdot=-2*rho[:,None,None]*p[:,0,None,None]*np.array([[A,C],[C,B]])
    transport_rhs=(kt[:,None]*current-div(spatial)
                   -np.einsum('nij,nj->ni',fdot,p[:,1:]))
    dp_transport=np.linalg.solve(kinetic,transport_rhs[...,None])[...,0]
    return dg,dk,-p,np.column_stack((dp0,dp_transport))

def add(state,step,scale): return tuple(x+scale*y for x,y in zip(state,step))
def rk4(state,h,dt):
    a=rhs(state,h); b=rhs(add(state,a,dt/2),h)
    c=rhs(add(state,b,dt/2),h); d=rhs(add(state,c,dt),h)
    return tuple(x+dt*(u+2*v+2*w+z)/6 for x,u,v,w,z in zip(state,a,b,c,d))

def initial(n):
    h=L/n; x=np.arange(n)*h
    f=np.column_stack((np.ones(n), .035*np.sin(2*np.pi*x/L),
                        .028*np.sin(2*np.pi*x/L+.4)))
    p=np.zeros_like(f)
    psi=np.ones(n)
    # Solve the periodic Lichnerowicz equation with constant isotropic K.
    for _ in range(500):
        g=np.repeat(psi[:,None]**4,3,axis=1)
        energy=matter(g,f,p,h)[0]
        ksq=24*np.pi*G*np.mean(psi**5*energy)/np.mean(psi**5)
        source=psi**5*(ksq/12-2*np.pi*G*energy)
        freq=2*np.pi*np.fft.fftfreq(n,d=h)
        shat=np.fft.fft(source); shat[0]=0
        candidate=np.real(np.fft.ifft(-shat/np.where(freq==0,1,freq**2)))+1
        if np.max(np.abs(candidate-psi))<1e-14: break
        psi=candidate
    g=np.repeat(psi[:,None]**4,3,axis=1)
    k=-np.sqrt(ksq)/3*g
    return x,(g,k,f,p)

def diagnostics(state,h):
    g,k,f,p=state; _,ham,mom,_=geometry(g,k,h)
    e,s,_,_=matter(g,f,p,h)
    ham-=16*np.pi*G*e; mom-=8*np.pi*G*s[:,0]
    return dict(hamiltonian_l2=float(np.sqrt(np.mean(ham**2))),
                momentum_l2=float(np.sqrt(np.mean(mom**2))),
                hamiltonian_max=float(np.max(np.abs(ham))),
                momentum_max=float(np.max(np.abs(mom))),
                rho_min=float(np.min(f[:,0])),metric_min=float(np.min(g)),
                kinetic_det_min=float(np.min((1+A*f[:,0]**2)*(1+B*f[:,0]**2)-(C*f[:,0]**2)**2)))

def run(n):
    x,state=initial(n); h=L/n
    start=diagnostics(state,h)
    steps=int(np.ceil(FINAL/(CFL*h))); dt=FINAL/steps
    for _ in range(steps): state=rk4(state,h,dt)
    end=diagnostics(state,h)
    return state,dict(n=n,dt=dt,steps=steps,initial=start,final=end)

def main():
    ns=[64,128,256]; solutions=[]; rows=[]
    for n in ns:
        state,row=run(n); rows.append(row); solutions.append(state)
    errors=[]
    for low,high in zip(solutions,solutions[1:]):
        errors.append(float(max(np.max(np.abs(u-v[::2])) for u,v in zip(low,high))))
    out=dict(scope='plane-symmetric diagonal 3+1 ADM, geodesic slicing, periodic domain; not toroidal BSSN',
             parameters=dict(G=G,L=L,t_end=FINAL,cfl=CFL,a=A,b=B,c=C),
             runs=rows,adjacent_resolution_max_errors=errors,
             observed_order=float(np.log2(errors[0]/errors[1])),
             gates=dict(positive_metric=all(r['final']['metric_min']>0 for r in rows),
                        positive_kinetic_matrix=all(r['final']['kinetic_det_min']>0 for r in rows),
                        finite_branch=all(r['final']['rho_min']>0 for r in rows),
                        refinement_error_decreases=errors[1]<errors[0],
                        final_hamiltonian_refines=all(rows[i+1]['final']['hamiltonian_l2'] < rows[i]['final']['hamiltonian_l2'] for i in range(2)),
                        final_momentum_refines=all(rows[i+1]['final']['momentum_l2'] < rows[i]['final']['momentum_l2'] for i in range(2))))
    path=Path(__file__).with_name('coupled_dee_adm_1d_checkpoint.json')
    path.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
