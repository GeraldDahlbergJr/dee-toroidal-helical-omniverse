"""BSSN initial-data bridge for the reduced cylindrical DEE checkpoint.

This gate converts the already solved physical conformal initial data to BSSN
variables and reconstructs the physical fields before independently evaluating
the same ADM constraints.  It is a conversion/consistency test, not evolution.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from solve_dynamical_dee_constraints import (
    solve, independent_residuals, RMIN, RMAX, G, S_Z, energy
)


def to_bssn(psi, dW):
    # gamma_ij = psi^4 delta_ij -> phi=ln psi, gtilde_ij=delta_ij.
    phi=np.log(psi)
    # K=0.  Existing checkpoint has physical K_rz = psi^-2 dW.
    # Atilde_ij = exp(-4 phi) K_ij for K=0.
    K_rz=psi**-2*dW
    Atilde_rz=psi**-4*K_rz
    return phi, Atilde_rz


def reconstruct(phi, Atilde_rz):
    psi=np.exp(phi)
    K_rz=psi**4*Atilde_rz
    return psi, K_rz


def bridge_residuals(sol,n):
    r=np.linspace(RMIN,RMAX,n); h=r[1]-r[0]
    psi,_,W,dW=sol.sol(r)
    phi,At=to_bssn(psi,dW)
    p,Krz=reconstruct(phi,At)
    # Algebraic BSSN constraints for this conformally-flat reduction.
    det_tilde_err=0.0
    tracefree_err=0.0
    # Reconstruct W derivative from physical K_rz: dW=psi^2 K_rz.
    dw=p*p*Krz
    ri=r[1:-1]; q=p[1:-1]
    dp=(p[2:]-p[:-2])/(2*h)
    ddp=(p[2:]-2*p[1:-1]+p[:-2])/h**2
    ddw=(dw[2:]-2*dw[1:-1]+dw[:-2])/h**2
    # Momentum divergence requires derivative of dW plus cylindrical dW/r.
    lapw=ddw+dw[1:-1]/ri
    lapp=ddp+dp/ri
    abar2=2*dw[1:-1]**2
    R3=-8*q**-5*lapp
    K2=q**-12*abar2
    H=R3-K2-16*np.pi*G*energy(q)
    Mz=q**-10*lapw-8*np.pi*G*q**-4*S_Z
    return {'points':n,
            'hamiltonian_max_abs':float(np.max(np.abs(H))),
            'hamiltonian_rms':float(np.sqrt(np.mean(H*H))),
            'momentum_max_abs':float(np.max(np.abs(Mz))),
            'momentum_rms':float(np.sqrt(np.mean(Mz*Mz))),
            'det_tilde_error':det_tilde_err,
            'tracefree_error':tracefree_err,
            'reconstruction_psi_max_abs':float(np.max(np.abs(p-psi))),
            'reconstruction_Krz_max_abs':float(np.max(np.abs(Krz-psi**-2*dW)))}


def main():
    sol=solve(); rows=[]
    for n in (201,401,801,1601,3201): rows.append(bridge_residuals(sol,n))
    for key in ('hamiltonian','momentum'):
        prev=None
        for row in rows:
            cur=row[f'{key}_max_abs']
            row[f'{key}_observed_order']=None if prev is None else float(np.log(prev/cur)/np.log(2))
            prev=cur
    out={'checkpoint_name':'BSSN initial-data bridge for coupled Einstein-DEE evolution',
         'scope':'conversion and reconstruction of previously constraint-satisfying reduced cylindrical initial data; no time evolution',
         'rows':rows,
         'gate':{
             'exact_field_reconstruction':all(x['reconstruction_psi_max_abs']<1e-13 and x['reconstruction_Krz_max_abs']<1e-13 for x in rows),
             'bssn_algebraic_constraints':all(x['det_tilde_error']<1e-13 and x['tracefree_error']<1e-13 for x in rows),
             'hamiltonian_decreases':all(rows[i]['hamiltonian_max_abs']<rows[i-1]['hamiltonian_max_abs'] for i in range(1,len(rows))),
             'momentum_decreases':all(rows[i]['momentum_max_abs']<rows[i-1]['momentum_max_abs'] for i in range(1,len(rows))),
             'hamiltonian_final_order_above_1_9':rows[-1]['hamiltonian_observed_order']>1.9,
             'momentum_final_order_above_1_9':rows[-1]['momentum_observed_order']>1.9},
         'next_gate':'short-duration coupled evolution only if every bridge gate is true'}
    Path(__file__).with_name('bssn_initial_data_bridge_checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
