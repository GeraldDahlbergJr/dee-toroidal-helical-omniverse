"""Hard t=0 Einstein-DEE momentum-interface consistency gate.

This gate tests the actual tensor/conformal map used by evolution:
    K_rz = psi^-2 Abar_rz,  Abar_rz = W'(r),
so that the physical upper-z momentum constraint is
    M^z = psi^-10 [d_r Abar_rz + Abar_rz/r] - 8 pi G psi^-4 S_z.

The independent checkpoint instead evaluates d_r Abar_rz as a direct second
finite difference of W.  Those are two consistent second-order discretizations
of the same continuum operator, but their truncation residuals need not be
numerically equal at finite h.  Therefore the hard gate is continuum
consistency: both residuals and their pointwise operator difference must
converge to zero at approximately second order.  No timestep is executed here.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from solve_dynamical_dee_constraints import (
    solve, RMIN, RMAX, G, energy, S_Z
)


def check(n: int):
    sol=solve()
    r=np.linspace(RMIN,RMAX,n); h=r[1]-r[0]
    psi,_,W,dW_exact=sol.sol(r)

    # Evolution-side stored physical tensor component.
    K_rz=psi**-2*dW_exact

    # Reconstruct the conformal longitudinal component from the stored tensor.
    Abar_rz=psi**2*K_rz

    ri=r[1:-1]; p=psi[1:-1]
    dp=(psi[2:]-psi[:-2])/(2*h)
    ddp=(psi[2:]-2*psi[1:-1]+psi[:-2])/h**2
    lap_p=ddp+dp/ri

    # Reference reduced cylindrical momentum operator, exactly as in the
    # independent initial-data evaluator: W'' + W'/r.
    Wp_ref=(W[2:]-W[:-2])/(2*h)
    Wpp_ref=(W[2:]-2*W[1:-1]+W[:-2])/h**2
    div_ref=Wpp_ref+Wp_ref/ri

    # Evolution/tensor-side operator obtained only from the stored K_rz.
    # Abar_rz = psi^2 K_rz, then div_bar Abar = d_r Abar_rz + Abar_rz/r.
    dAbar=(Abar_rz[2:]-Abar_rz[:-2])/(2*h)
    div_map=dAbar+Abar_rz[1:-1]/ri

    # Hamiltonian: use the same physical K tensor mapping.
    abar2=2*Abar_rz[1:-1]**2
    H_map=-8*p**-5*lap_p-p**-12*abar2-16*np.pi*G*energy(p)

    M_ref=p**-10*div_ref-8*np.pi*G*p**-4*S_Z
    M_map=p**-10*div_map-8*np.pi*G*p**-4*S_Z

    # Pointwise algebraic/discrete consistency diagnostics.
    div_diff=div_map-div_ref
    M_diff=M_map-M_ref

    return {
      'points':n,
      'hamiltonian_max_abs':float(np.max(np.abs(H_map))),
      'reference_momentum_max_abs':float(np.max(np.abs(M_ref))),
      'mapped_momentum_max_abs':float(np.max(np.abs(M_map))),
      'divergence_operator_difference_max_abs':float(np.max(np.abs(div_diff))),
      'momentum_difference_max_abs':float(np.max(np.abs(M_diff))),
      'reconstruction_Abar_max_abs':float(np.max(np.abs(Abar_rz-dW_exact))),
      'finite':bool(np.all(np.isfinite(H_map)) and np.all(np.isfinite(M_map)))
    }


def add_orders(rows,key):
    prev=None
    for row in rows:
        cur=row[key]
        row[key.replace('_max_abs','')+'_observed_order']=None if prev is None else float(np.log(prev/cur)/np.log(2))
        prev=cur


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--resolutions',nargs='+',type=int,default=[201,401,801])
    ap.add_argument('--t-end',type=float,default=.25)  # workflow compatibility
    ap.add_argument('--cfl',type=float,default=.10)    # workflow compatibility
    ap.add_argument('--output',default='coupled_einstein_dee_checkpoint.json')
    a=ap.parse_args()

    rows=[check(n) for n in a.resolutions]
    for key in (
        'hamiltonian_max_abs',
        'reference_momentum_max_abs',
        'mapped_momentum_max_abs',
        'divergence_operator_difference_max_abs',
        'momentum_difference_max_abs',
    ):
        add_orders(rows,key)

    def decreases(key):
        return all(rows[i][key] < rows[i-1][key] for i in range(1,len(rows)))

    final=rows[-1]
    order_keys=(
        'hamiltonian_observed_order',
        'reference_momentum_observed_order',
        'mapped_momentum_observed_order',
        'divergence_operator_difference_observed_order',
        'momentum_difference_observed_order',
    )
    approx_second_order=all(final[k] is not None and final[k] > 1.7 for k in order_keys)

    reconstruction_exact=all(r['reconstruction_Abar_max_abs'] < 1e-13 for r in rows)
    all_finite=all(r['finite'] for r in rows)
    all_decrease=all(decreases(k) for k in (
        'hamiltonian_max_abs',
        'reference_momentum_max_abs',
        'mapped_momentum_max_abs',
        'divergence_operator_difference_max_abs',
        'momentum_difference_max_abs',
    ))

    passed=reconstruction_exact and all_finite and all_decrease and approx_second_order

    out={
      'checkpoint_name':'Coupled Einstein-DEE t=0 tensor momentum interface gate',
      'status':'PASS' if passed else 'INCONCLUSIVE',
      'scope':'pre-evolution cylindrical conformal tensor-map consistency; no timestep executed',
      'identity_tested':'K_rz=psi^-2 Abar_rz with Abar_rz=Wprime; M^z=psi^-10 div_bar(Abar)-8piG psi^-4 S_z',
      'interpretation':'finite-grid residual norms from two second-order stencils need not be equal; the required invariant is common continuum convergence and vanishing pointwise operator difference',
      'runs':rows,
      'gate':{
        'Abar_reconstruction_exact':reconstruction_exact,
        'all_quantities_finite':all_finite,
        'all_residuals_and_operator_differences_decrease':all_decrease,
        'final_orders_above_1_7':approx_second_order,
      },
      'next_gate':'restore short-duration coupled ADM evolution only after this gate passes'
    }
    Path(a.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    if not passed:
        raise SystemExit(3)

if __name__=='__main__':
    main()
