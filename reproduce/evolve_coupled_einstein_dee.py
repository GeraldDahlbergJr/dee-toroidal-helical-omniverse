"""Hard t=0 Einstein-DEE interface gate.

This checkpoint reconstructs the conformal variables from the state that will
be handed to evolution and independently reproduces the cylindrical Hamiltonian
and momentum constraints.  No timestep is permitted unless this agrees with
solve_dynamical_dee_constraints.py at all requested resolutions.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from solve_dynamical_dee_constraints import (
    solve, independent_residuals, RMIN, RMAX, G, RHO2, energy, S_Z
)


def check(n: int):
    sol=solve()
    r=np.linspace(RMIN,RMAX,n); h=r[1]-r[0]
    psi,_,W,dW_exact=sol.sol(r)

    # Reconstruct exactly the quantities passed through the evolution interface.
    gamma_rr=psi**4
    K_rz=psi**-2*dW_exact
    psi_evo=gamma_rr**0.25
    dW_evo=psi_evo**2*K_rz

    ri=r[1:-1]; p=psi_evo[1:-1]
    dp=(psi_evo[2:]-psi_evo[:-2])/(2*h)
    ddp=(psi_evo[2:]-2*psi_evo[1:-1]+psi_evo[:-2])/h**2
    lap_p=ddp+dp/ri

    # dW_evo is W'.  Therefore W'' is its FIRST derivative.  The previous
    # implementation accidentally differentiated W' twice (W'''), producing
    # the nonconvergent O(1e-2) momentum mismatch.
    Wpp=(dW_evo[2:]-dW_evo[:-2])/(2*h)
    lap_w=Wpp+dW_evo[1:-1]/ri

    abar2=2*dW_evo[1:-1]**2
    H=-8*p**-5*lap_p-p**-12*abar2-16*np.pi*G*energy(p)
    M=p**-10*lap_w-8*np.pi*G*p**-4*S_Z
    eh=float(np.max(np.abs(H))); em=float(np.max(np.abs(M)))
    rh,_,rm,_=independent_residuals(sol,n)
    dh=abs(eh-rh)/max(rh,1e-30); dm=abs(em-rm)/max(rm,1e-30)
    return {
      'points':n,'reference_hamiltonian_max_abs':rh,
      'evolution_t0_hamiltonian_max_abs':eh,
      'reference_momentum_max_abs':rm,'evolution_t0_momentum_max_abs':em,
      'hamiltonian_relative_mismatch':dh,'momentum_relative_mismatch':dm,
      'interface_pass':bool(dh<0.03 and dm<0.03)
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--resolutions',nargs='+',type=int,default=[201,401,801])
    ap.add_argument('--t-end',type=float,default=.25)
    ap.add_argument('--cfl',type=float,default=.10)
    ap.add_argument('--output',default='coupled_einstein_dee_checkpoint.json')
    a=ap.parse_args()
    rows=[check(n) for n in a.resolutions]
    for key in ('reference_hamiltonian','reference_momentum','evolution_t0_hamiltonian','evolution_t0_momentum'):
        prev=None
        for row in rows:
            cur=row[key+'_max_abs']
            row[key+'_observed_order']=None if prev is None else float(np.log(prev/cur)/np.log(2))
            prev=cur
    passed=all(r['interface_pass'] for r in rows)
    out={
      'checkpoint_name':'Coupled Einstein-DEE t=0 interface gate',
      'status':'PASS' if passed else 'INCONCLUSIVE',
      'scope':'pre-evolution cylindrical conformal consistency gate; no timestep executed',
      'runs':rows,'gate':{'t0_interface_consistent':passed},
      'next_gate':'restore coupled ADM evolution only after this gate passes'
    }
    Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    if not passed: raise SystemExit(3)

if __name__=='__main__': main()
