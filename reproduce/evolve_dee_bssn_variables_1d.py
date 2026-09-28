"""Plane-symmetric BSSN evolution with an evolved conformal connection.

Conformal variables (phi, tilde_gamma_ii, K, tilde_A_ii, tilde_Gamma^x)
are independently evolved. The conformal Ricci operator uses the evolved
connection. This symmetry-restricted implementation is not a full
multidimensional or toroidal BSSN implementation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from evolve_dee_adm_1d import (G,L,FINAL,CFL,A,B,C,dx,geometry,matter,rhs as adm_rhs,
                               initial as adm_initial,diagnostics as adm_diagnostics)


def to_bssn(state,h):
    g,k,f,p=state
    phi=np.log(np.prod(g,axis=1))/12
    scale=np.exp(4*phi)
    tg=g/scale[:,None]
    kt=np.sum(k/g,axis=1)
    ta=(k-g*kt[:,None]/3)/scale[:,None]
    contracted=dx(tg[:,0],h)/tg[:,0]**2
    return phi,tg,kt,ta,contracted,f,p


def to_adm(state):
    phi,tg,kt,ta,_,f,p=state
    g=np.exp(4*phi)[:,None]*tg
    k=np.exp(4*phi)[:,None]*(ta+tg*kt[:,None]/3)
    return g,k,f,p


def conformal_diagnostics(state,h):
    phi,tg,kt,ta,conn,f,p=state
    adm=to_adm(state)
    d=adm_diagnostics(adm,h)
    determinant=np.prod(tg,axis=1)
    atrace=np.sum(ta/tg,axis=1)
    connection=dx(tg[:,0],h)/tg[:,0]**2
    d.update(conformal_det_max=float(np.max(np.abs(determinant-1))),
             conformal_trace_max=float(np.max(np.abs(atrace))),
             connection_constraint_l2=float(np.sqrt(np.mean((conn-connection)**2))),
             connection_constraint_max=float(np.max(np.abs(conn-connection))))
    return d


def conformal_ricci(phi,tg,conn,h):
    """Physical diagonal R_ij = tilde R_ij + R^phi_ij in 1-D BSSN.

    tilde R is obtained from the conformal metric and its contracted
    connection, including the evolved connection in both its derivative
    and its undifferentiated Christoffel term. With the connection constraint
    satisfied, this reduces to the conformal geometric Ricci tensor.
    """
    tr,_,_,_=geometry(tg,np.zeros_like(tg),h)
    td=dx(tg,h)
    ch=np.column_stack((td[:,0]/(2*tg[:,0]),
                        -td[:,1]/(2*tg[:,0]),
                        -td[:,2]/(2*tg[:,0])))
    derived=td[:,0]/tg[:,0]**2
    delta=conn-derived
    tr+=tg[:,0,None]*ch*delta[:,None]
    tr[:,0]+=tg[:,0]*dx(delta,h)
    dp=dx(phi,h); ddp=(np.roll(phi,-1)-2*phi+np.roll(phi,1))/h**2
    cov=np.column_stack((ddp-ch[:,0]*dp,-ch[:,1]*dp,-ch[:,2]*dp))
    lap=np.sum(cov/tg,axis=1)
    grad2=dp**2/tg[:,0]
    rphi=-2*cov-2*tg*lap[:,None]-4*tg*grad2[:,None]
    rphi[:,0]+=4*dp**2
    return tr+rphi


def bssn_rhs(state,h):
    phi,tg,kt,ta,conn,f,p=state
    g,k,_,_=to_adm(state)
    if np.min(g)<=0 or np.min(tg)<=0: raise FloatingPointError('invalid metric signature')
    ric=conformal_ricci(phi,tg,conn,h)
    e,s,stress,trace=matter(g,f,p,h)
    sd=np.diagonal(stress,axis1=1,axis2=2)
    scale=np.exp(4*phi)
    aa=np.sum((ta/tg)**2,axis=1)
    # Hamiltonian-constraint substituted BSSN trace equation.
    dkt=aa+kt**2/3+4*np.pi*G*(e+trace)
    q=ric-8*np.pi*G*sd
    qtrace=np.sum(q/g,axis=1)
    dta=(q-g*qtrace[:,None]/3)/scale[:,None]
    dta+=kt[:,None]*ta-2*ta*ta/tg
    # Momentum-constraint substituted conformal connection equation.
    gder=dx(tg,h)
    chx=gder[:,0]/(2*tg[:,0]); chy=-gder[:,1]/(2*tg[:,0]); chz=-gder[:,2]/(2*tg[:,0])
    christoffel_a=(chx*ta[:,0]/tg[:,0]**2
                   +chy*ta[:,1]/tg[:,1]**2+chz*ta[:,2]/tg[:,2]**2)
    dconn=(2*christoffel_a+12*(ta[:,0]/tg[:,0]**2)*dx(phi,h)
           -4/(3*tg[:,0])*dx(kt,h)-16*np.pi*G*s[:,0]/tg[:,0])
    # Reuse the identical nonlinear DEE equations to isolate geometric differences.
    _,_,df,dp=adm_rhs((g,k,f,p),h)
    return -kt/6,-2*ta,dkt,dta,dconn,df,dp


def add(state,step,dt): return tuple(x+dt*y for x,y in zip(state,step))
def rk4(state,h,dt):
    a=bssn_rhs(state,h); b=bssn_rhs(add(state,a,dt/2),h)
    c=bssn_rhs(add(state,b,dt/2),h); d=bssn_rhs(add(state,c,dt),h)
    return tuple(x+dt*(u+2*v+2*w+z)/6 for x,u,v,w,z in zip(state,a,b,c,d))


def run(n):
    h=L/n; _,adm=adm_initial(n)
    state=to_bssn(adm,h); initial=conformal_diagnostics(state,h)
    steps=int(np.ceil(FINAL/(CFL*h))); dt=FINAL/steps
    for _ in range(steps): state=rk4(state,h,dt)
    return state,dict(n=n,dt=dt,steps=steps,initial=initial,
                      final=conformal_diagnostics(state,h))


def main():
    states=[]; rows=[]; comparisons=[]
    for n in (64,128,256):
        state,row=run(n); states.append(state); rows.append(row)
        # Full ADM evolution is an independent geometric-RHS comparison.
        from evolve_dee_adm_1d import run as adm_run
        other,_=adm_run(n)
        comparisons.append(float(max(np.max(np.abs(x-y)) for x,y in
                                     zip(to_adm(state),other))))
    diffs=[]
    for low,high in zip(states,states[1:]):
        diffs.append(float(max(np.max(np.abs(x-y[::2])) for x,y in
                               zip(to_adm(low),to_adm(high)))))
    out=dict(scope='plane-symmetric diagonal BSSN conformal variables, evolved connection and conformal Ricci; not full 3-D or toroidal geometry',
             parameters=dict(G=G,L=L,t_end=FINAL,cfl=CFL),runs=rows,
             adjacent_resolution_max_errors=diffs,
             observed_order=float(np.log2(diffs[0]/diffs[1])),
             same_resolution_adm_max_differences=comparisons,
             gates=dict(hamiltonian_refines=all(rows[i+1]['final']['hamiltonian_l2']<rows[i]['final']['hamiltonian_l2'] for i in range(2)),
                        momentum_refines=all(rows[i+1]['final']['momentum_l2']<rows[i]['final']['momentum_l2'] for i in range(2)),
                        connection_refines=all(rows[i+1]['final']['connection_constraint_l2']<rows[i]['final']['connection_constraint_l2'] for i in range(2)),
                        positive_metric=all(row['final']['metric_min']>0 for row in rows)))
    Path(__file__).with_name('coupled_dee_bssn_variables_1d_checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
