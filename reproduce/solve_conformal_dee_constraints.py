"""Solve a first nonlinear DEE 3+1 constraint gate.

Restricted gate:
* local cylindrical annulus r in [r_min,R]
* conformally flat spatial metric gamma_ij = psi^4 delta_ij
* maximal/time-symmetric K_ij = 0
* Pi_rho=Pi_Theta=Pi_Psi=0, hence S_i=0 and M_i=0 identically
* rho=sqrt(0.945), Theta=m_Theta phi+k_Theta z,
  Psi=m_Psi phi+k_Psi z

The Hamiltonian constraint is solved nonlinearly for psi.  This is genuine
constraint-satisfying initial data for this restricted ansatz, not a full
3-D evolution and not a proof of nonlinear DEE stability.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp

G=0.005
RHO2=0.945
RHO=np.sqrt(RHO2)
LAMBDA_V=4.0; V0=1.0; A=0.30; B=0.20; C=0.10
RMIN=0.5; RMAX=8.0
M_THETA=1.0; M_PSI=-1.0; K_THETA=0.25; K_PSI=0.20

KT=1.0+A*RHO2; KP=1.0+B*RHO2; MIX=C*RHO2
POT=0.25*LAMBDA_V*(RHO2-V0*V0)**2

def q(r):
    xt=K_THETA**2+M_THETA**2/r**2
    xp=K_PSI**2+M_PSI**2/r**2
    cross=K_THETA*K_PSI+M_THETA*M_PSI/r**2
    return KT*xt+KP*xp+2.0*MIX*cross

def energy(r,psi):
    return POT+0.5*psi**-4*q(r)

def ode(r,y):
    psi,dpsi=y
    # Delta psi = -2 pi G psi^5 E; cylindrical Delta=d2+dr/r.
    rhs=-2.0*np.pi*G*(POT*psi**5+0.5*q(r)*psi)
    return np.vstack((dpsi,rhs-dpsi/r))

def bc(ya,yb):
    return np.array([ya[1],yb[0]-1.0])

def solve():
    r=np.linspace(RMIN,RMAX,400)
    y=np.vstack((np.ones_like(r),np.zeros_like(r)))
    sol=solve_bvp(ode,bc,r,y,tol=1e-11,max_nodes=30000)
    if sol.status != 0:
        raise RuntimeError(sol.message)
    return sol

def residuals(sol,n):
    r=np.linspace(RMIN,RMAX,n); psi=sol.sol(r)[0]; h=r[1]-r[0]
    d=(psi[2:]-psi[:-2])/(2*h)
    dd=(psi[2:]-2*psi[1:-1]+psi[:-2])/h**2
    lap=dd+d/r[1:-1]
    e=energy(r[1:-1],psi[1:-1])
    R3=-8.0*psi[1:-1]**-5*lap
    H=R3-16.0*np.pi*G*e
    # K_ij=0 and all canonical momenta vanish, so S_i=0 exactly.
    M=np.zeros((n-2,3))
    return float(np.max(np.abs(H))),float(np.sqrt(np.mean(H*H))),float(np.max(np.abs(M)))

def main():
    sol=solve(); rows=[]; prev=None
    for n in (201,401,801,1601,3201):
        hmax,hrms,mmax=residuals(sol,n)
        order=None if prev is None else float(np.log(prev/hmax)/np.log(2.0))
        rows.append({'points':n,'hamiltonian_max_abs':hmax,'hamiltonian_rms':hrms,
                     'hamiltonian_observed_order':order,'momentum_max_abs':mmax})
        prev=hmax
    out={
      'scope':'restricted local cylindrical-annulus conformally-flat time-symmetric nonlinear 3+1 DEE initial-data gate',
      'not_claimed':'full 3-D Einstein-matter evolution or nonlinear stability',
      'parameters':{'G':G,'rho_squared':RHO2,'r_min':RMIN,'r_max':RMAX,
                    'm_theta':M_THETA,'m_psi':M_PSI,'k_theta':K_THETA,'k_psi':K_PSI},
      'solution':{'psi_min':float(np.min(sol.y[0])),'psi_max':float(np.max(sol.y[0])),
                  'inner_neumann_abs':float(abs(sol.y[1,0])),'outer_dirichlet_abs':float(abs(sol.y[0,-1]-1.0))},
      'refinement':rows,
      'gate':{'hamiltonian_residual_decreases':all(rows[i]['hamiltonian_max_abs']<rows[i-1]['hamiltonian_max_abs'] for i in range(1,len(rows))),
              'asymptotic_order_above_1_9':rows[-1]['hamiltonian_observed_order']>1.9,
              'momentum_constraint_exact_in_restricted_time_symmetric_sector':all(x['momentum_max_abs']==0.0 for x in rows)}
    }
    path=Path(__file__).with_name('conformal_dee_constraint_checkpoint.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
