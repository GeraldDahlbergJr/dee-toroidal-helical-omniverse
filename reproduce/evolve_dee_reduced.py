"""Reduced nonlinear DEE matter-evolution gate.

This is deliberately *not* advertised as a full BSSN/GH 3+1 Einstein-matter
code.  It evolves the action-derived DEE radial field on a fixed local
background as a first dynamical matter/recovery checkpoint.  Full metric
backreaction remains a subsequent gate.

Equation tested (restricted homogeneous transport invariants):
    rho_tt = rho_rr + rho_r/r - dVeff/drho
with
    dVeff/drho = lambda*rho*(rho^2-v^2)
                  + rho*(a X_theta+b X_psi+2 c X_cross).
The equilibrium rho_* is therefore the archived finite DEE branch.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

LAM=4.0; V=1.0; A=.30; B=.20; C=.10
XT=.60; XP=.40; XC=-.20
Q=A*XT+B*XP+2*C*XC
RHO2=V*V-Q/LAM
RHO=np.sqrt(RHO2)
RMIN=.5; RMAX=8.0
TEND=8.0; CFL=.20
AMP=.02; CENTER=3.0; WIDTH=.55

def force(rho):
    return LAM*rho*(rho*rho-V*V)+Q*rho

def rhs(r,rho):
    h=r[1]-r[0]
    out=np.zeros_like(rho)
    d=(rho[2:]-rho[:-2])/(2*h)
    dd=(rho[2:]-2*rho[1:-1]+rho[:-2])/h**2
    out[1:-1]=dd+d/r[1:-1]-force(rho[1:-1])
    # homogeneous Neumann boundaries for this local annulus gate
    out[0]=2*(rho[1]-rho[0])/h**2-force(rho[0])
    out[-1]=2*(rho[-2]-rho[-1])/h**2-force(rho[-1])
    return out

def run(n):
    r=np.linspace(RMIN,RMAX,n); h=r[1]-r[0]; dt=CFL*h
    steps=int(np.ceil(TEND/dt)); dt=TEND/steps
    rho=RHO+AMP*np.exp(-((r-CENTER)/WIDTH)**2)
    pi=np.zeros_like(rho)
    # velocity Verlet / leapfrog, second order
    acc=rhs(r,rho); pi += .5*dt*acc
    maxdev0=float(np.max(np.abs(rho-RHO)))
    maxdev=maxdev0
    for k in range(steps):
        rho += dt*pi
        acc=rhs(r,rho)
        if k == steps-1:
            pi += .5*dt*acc
        else:
            pi += dt*acc
        maxdev=max(maxdev,float(np.max(np.abs(rho-RHO))))
    return r,rho,pi,{'points':n,'dt':dt,'steps':steps,
                     'initial_max_deviation':maxdev0,
                     'final_max_deviation':float(np.max(np.abs(rho-RHO))),
                     'max_deviation_over_run':maxdev,
                     'rho_min':float(np.min(rho)),'rho_max':float(np.max(rho))}

def main():
    ns=(201,401,801,1601)
    sols=[]; rows=[]
    for n in ns:
        r,rho,pi,row=run(n); sols.append((r,rho)); rows.append(row)
    # compare coarse final solutions to finest sampled at identical nodes
    rf,uf=sols[-1]
    errs=[]
    for r,u in sols[:-1]:
        ref=np.interp(r,rf,uf)
        errs.append(float(np.max(np.abs(u-ref))))
    orders=[]
    for i in range(len(errs)-1):
        orders.append(float(np.log(errs[i]/errs[i+1])/np.log(2)))
    out={'scope':'reduced nonlinear DEE matter evolution on fixed local background; not full Einstein-matter evolution',
         'rho_star_squared':RHO2,'rho_star':RHO,'perturbation_amplitude':AMP,'t_end':TEND,
         'runs':rows,'coarse_to_finest_max_errors':errs,'observed_orders':orders,
         'gate':{'finite_branch_remains_positive':all(x['rho_min']>0 for x in rows),
                 'bounded_over_interval':all(x['max_deviation_over_run']<0.10 for x in rows),
                 'refinement_error_decreases':all(errs[i+1]<errs[i] for i in range(len(errs)-1))}}
    p=Path(__file__).with_name('reduced_dee_evolution_checkpoint.json')
    p.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
