"""Coupled nonlinear DEE initial-data constraint checkpoint.

Reduced local cylindrical-annulus gate.  Unlike the preceding time-symmetric
checkpoint this uses nonzero transport-field normal momenta, S_z != 0, and a
nonzero conformal vector potential W_z, hence K_rz != 0.

Conformal ansatz: gamma_ij=psi^4 delta_ij, K=0,
K^ij=psi^-10 (L W)^ij.  The active vector component is W_z(r).
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp

G=0.005; RHO2=0.945
KT=1.0+0.30*RHO2; KP=1.0+0.20*RHO2; MIX=0.10*RHO2
POT=0.25*4.0*(RHO2-1.0)**2
RMIN=0.5; RMAX=8.0
K_THETA=0.25; K_PSI=0.20
PI_THETA=0.12; PI_PSI=-0.04
Q=KT*K_THETA**2+KP*K_PSI**2+2*MIX*K_THETA*K_PSI
PI_KIN=KT*PI_THETA**2+KP*PI_PSI**2+2*MIX*PI_THETA*PI_PSI
S_Z=(KT*PI_THETA*K_THETA+KP*PI_PSI*K_PSI
     +MIX*(PI_THETA*K_PSI+PI_PSI*K_THETA))

def energy(psi):
    return POT+0.5*PI_KIN+0.5*psi**-4*Q

def ode(r,y):
    psi,dpsi,W,dW=y
    abar2=2.0*dW*dW
    ddpsi=-(1/8)*psi**-7*abar2-2*np.pi*G*psi**5*energy(psi)-dpsi/r
    ddW=8*np.pi*G*psi**6*S_Z-dW/r
    return np.vstack((dpsi,ddpsi,dW,ddW))

def bc(ya,yb):
    return np.array([ya[1],yb[0]-1.0,ya[3],yb[2]])

def solve():
    r=np.linspace(RMIN,RMAX,500)
    y=np.vstack((np.ones_like(r),np.zeros_like(r),np.zeros_like(r),np.zeros_like(r)))
    sol=solve_bvp(ode,bc,r,y,tol=1e-11,max_nodes=50000)
    if sol.status != 0: raise RuntimeError(sol.message)
    return sol

def independent_residuals(sol,n):
    r=np.linspace(RMIN,RMAX,n); h=r[1]-r[0]
    psi,_,W,_=sol.sol(r); ri=r[1:-1]; p=psi[1:-1]
    dp=(psi[2:]-psi[:-2])/(2*h)
    ddp=(psi[2:]-2*psi[1:-1]+psi[:-2])/h**2
    dw=(W[2:]-W[:-2])/(2*h)
    ddw=(W[2:]-2*W[1:-1]+W[:-2])/h**2
    lap_p=ddp+dp/ri; lap_w=ddw+dw/ri
    abar2=2*dw*dw
    R3=-8*p**-5*lap_p
    K2=p**-12*abar2
    H=R3-K2-16*np.pi*G*energy(p)
    # Physical upper-z momentum residual: psi^-10 div_bar(Abar)-8piG psi^-4 S_z.
    Mz=p**-10*lap_w-8*np.pi*G*p**-4*S_Z
    return (float(np.max(np.abs(H))),float(np.sqrt(np.mean(H*H))),
            float(np.max(np.abs(Mz))),float(np.sqrt(np.mean(Mz*Mz))))

def main():
    sol=solve(); rows=[]; prev_h=prev_m=None
    for n in (201,401,801,1601,3201):
        hm,hr,mm,mr=independent_residuals(sol,n)
        ho=None if prev_h is None else float(np.log(prev_h/hm)/np.log(2))
        mo=None if prev_m is None else float(np.log(prev_m/mm)/np.log(2))
        rows.append({'points':n,'hamiltonian_max_abs':hm,'hamiltonian_rms':hr,
                     'hamiltonian_observed_order':ho,'momentum_max_abs':mm,
                     'momentum_rms':mr,'momentum_observed_order':mo})
        prev_h,prev_m=hm,mm
    r=np.linspace(RMIN,RMAX,5000); psi,_,_,dW=sol.sol(r)
    out={'scope':'reduced local cylindrical-annulus conformal dynamical nonlinear DEE initial-data gate',
         'not_claimed':'full 3-D evolution, global toroidal solution, or nonlinear stability',
         'parameters':{'G':G,'rho_squared':RHO2,'pi_theta':PI_THETA,'pi_psi':PI_PSI,
                       'k_theta':K_THETA,'k_psi':K_PSI,'S_z':S_Z},
         'solution':{'psi_min':float(np.min(psi)),'psi_max':float(np.max(psi)),
                     'max_abs_S_i':float(abs(S_Z)),
                     'max_abs_K_rz':float(np.max(np.abs(psi**-2*dW)))},
         'refinement':rows,
         'gate':{'matter_momentum_nonzero':abs(S_Z)>0,
                 'extrinsic_curvature_nonzero':float(np.max(np.abs(psi**-2*dW)))>0,
                 'hamiltonian_residual_decreases':all(rows[i]['hamiltonian_max_abs']<rows[i-1]['hamiltonian_max_abs'] for i in range(1,len(rows))),
                 'momentum_residual_decreases':all(rows[i]['momentum_max_abs']<rows[i-1]['momentum_max_abs'] for i in range(1,len(rows))),
                 'hamiltonian_asymptotic_order_above_1_9':rows[-1]['hamiltonian_observed_order']>1.9,
                 'momentum_asymptotic_order_above_1_9':rows[-1]['momentum_observed_order']>1.9}}
    Path(__file__).with_name('dynamical_dee_constraint_checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
