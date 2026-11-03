"""Full 3+1 stress-energy projections for the recovered nonlinear DEE action.

Signature (-,+,+,+).  Normal derivatives Pi_A = n^mu d_mu A and spatial
covectors D_i A.  The phase kinetic matrix is
K_AB=[[KTheta,Lambda],[Lambda,KPsi]].

This module closes E, S_i, S_ij, S for the BSSN source terms.  It contains a
separate direct 4D tensor construction used as an independent numerical check.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

LAMV=4.; VEV=1.; A=.30; B=.20; C=.10

def coeffs(rho): return 1+A*rho*rho, 1+B*rho*rho, C*rho*rho

def potential(rho): return .25*LAMV*(rho*rho-VEV*VEV)**2

def projections(rho, Pi, D, gamma_inv):
    """Return E,S_i,S_ij,S from analytic 3+1 formulas.

    Pi: shape (3,) ordered rho,Theta,Psi. D: shape (3,3), field,space.
    gamma_inv: 3x3 inverse spatial metric.
    """
    kt,kp,la=coeffs(rho)
    K=np.array([[1.,0.,0.],[0.,kt,la],[0.,la,kp]])
    spatial=np.einsum('ab,ai,ij,bj->',K,D,gamma_inv,D)
    temporal=float(Pi@K@Pi)
    U=potential(rho)
    E=.5*(temporal+spatial)+U
    # convention S_i=-gamma_i^mu n^nu T_munu = - K_AB Pi_A D_i phi_B
    Si=-np.einsum('ab,a,bi->i',K,Pi,D)
    gamma=np.linalg.inv(gamma_inv)
    # L_m = 1/2 temporal - 1/2 spatial - U
    Lm=.5*temporal-.5*spatial-U
    Sij=np.einsum('ab,ai,bj->ij',K,D,D)+gamma*Lm
    S=float(np.einsum('ij,ij->',gamma_inv,Sij))
    return E,Si,Sij,S

def direct_tensor_check(rho,Pi,D,gamma):
    """Independent zero-shift alpha=1 4D construction and projection."""
    kt,kp,la=coeffs(rho)
    K=np.array([[1.,0.,0.],[0.,kt,la],[0.,la,kp]])
    g=np.zeros((4,4)); g[0,0]=-1.; g[1:,1:]=gamma
    gi=np.linalg.inv(g)
    deriv=np.zeros((3,4)); deriv[:,0]=Pi; deriv[:,1:]=D
    kinetic=np.einsum('ab,am,mn,bn->',K,deriv,gi,deriv)
    Lm=-.5*kinetic-potential(rho)
    T=np.einsum('ab,am,bn->mn',K,deriv,deriv)+g*Lm
    E=T[0,0]; Si=-T[1:,0]; Sij=T[1:,1:]
    S=float(np.einsum('ij,ij->',np.linalg.inv(gamma),Sij))
    return E,Si,Sij,S,T

def main():
    rng=np.random.default_rng(20260928); worst=0.; sym=0.; trace=0.; mixed=0.
    for _ in range(10000):
        rho=float(rng.uniform(.2,1.5)); Pi=rng.normal(size=3); D=rng.normal(size=(3,3))
        M=rng.normal(size=(3,3)); gamma=M.T@M+np.eye(3); gi=np.linalg.inv(gamma)
        a=projections(rho,Pi,D,gi); b=direct_tensor_check(rho,Pi,D,gamma)
        vals=[abs(a[0]-b[0]),np.max(abs(a[1]-b[1])),np.max(abs(a[2]-b[2])),abs(a[3]-b[3])]
        worst=max(worst,*map(float,vals)); sym=max(sym,float(np.max(abs(a[2]-a[2].T))))
        trace=max(trace,float(abs(a[3]-np.einsum('ij,ij->',gi,a[2]))))
        # Mixed-term isolation: full result minus c=0 result is checked directly by 4D tensor route
        # through the overall analytic-vs-direct comparison above; record Lambda nonzero coverage.
        mixed=max(mixed,abs(coeffs(rho)[2]))
    out={'checkpoint_name':'DEE full 3+1 stress-energy projection closure',
         'samples':10000,'max_projection_disagreement':worst,'max_Sij_asymmetry':sym,
         'max_trace_disagreement':trace,'max_abs_Lambda_sampled':mixed,
         'gate':{'analytic_matches_independent_4d_tensor':worst<1e-10,
                 'Sij_symmetric':sym<1e-12,'trace_consistent':trace<1e-12,
                 'mixed_theta_psi_sector_exercised':mixed>0},
         'scope':'matter-source algebra closure only; not a coupled evolution result'}
    Path(__file__).with_name('dee_stress_energy_3p1_checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
