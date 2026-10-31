#!/usr/bin/env python3
"""Experimental C3 alpha/beta helical coupling with fresh nonlinear ADM solves.

No Run43 frozen field is changed. This is a modified action-derived matter
source model, NOT the frozen Run43 formulation or its 12-check validation.
"""
import hashlib
import json
import math
from pathlib import Path
import numpy as np

from dee_stress_energy_3p1 import coeffs
from validate_toroidal_helical_dee_source import RMIN, RMAX, M, RHO, PI, KTH, KPS
from solve_toroidal_helical_coupled_constraints import Grid
from diagnose_toroidal_boundary_convergence import Audit

HERE = Path(__file__).resolve().parent
EPS = 0.002
OMEGA = 0.4
WIDTH = 0.06
LEVELS = (8, 16, 32)


def ramp(r):
    """Septic C3 switch: beta at inner edge, two-quadrature alpha in bulk."""
    u = np.clip((r - RMIN) / WIDTH, 0., 1.)
    active = (u > 0.) & (u < 1.)
    s = 35*u**4 - 84*u**5 + 70*u**6 - 20*u**7
    ds = np.where(active, (140*u**3 - 420*u**4 + 420*u**5 - 140*u**6)/WIDTH, 0.)
    return s, ds


def matter_arrays(r,t,p,q,basis,shape):
    er,et,ep = basis
    s,ds=ramp(r)
    chi=t-M*p
    direction=et/r-M*ep/q
    kinetic=np.array([[1.,0.,0.],[0.,coeffs(RHO)[0],coeffs(RHO)[2]],
                       [0.,coeffs(RHO)[2],coeffs(RHO)[1]]])
    D=np.zeros((3,3)+shape)
    D[1]=(KTH-EPS*np.sin(chi))*direction
    D[2]=(KPS+EPS*s*np.cos(chi))*direction+EPS*ds*np.sin(chi)*er
    velocities=np.zeros((3,)+shape)
    for i in range(3):
        velocities[i]=PI[i]
    velocities[1]-=EPS*OMEGA*np.sin(chi)
    velocities[2]+=EPS*OMEGA*s*np.cos(chi)
    temporal=np.einsum('ab,a...,b...->...',kinetic,velocities,velocities)
    spatial=np.einsum('ab,ai...,bi...->...',kinetic,D,D)
    momentum=-np.einsum('ab,a...,bi...->i...',kinetic,velocities,D)
    return temporal,spatial,momentum


class PhaseGrid(Grid):
    def __init__(self,nr,nt,nph):
        super().__init__(nr,nt,nph)
        temp,spatial,momentum=matter_arrays(
            self.r,self.t,self.p,self.q,self.basis,(self.size,))
        self.temporal=temp
        self.spatial=spatial
        self.S=momentum


class PhaseAudit(Audit):
    def __init__(self,shape):
        super().__init__(shape)
        er,et,ep=self.er,self.et,self.ep
        basis=(er,et,ep)
        r=np.broadcast_to(self.r,shape)
        t=np.broadcast_to(self.t,shape)
        p=np.broadcast_to(self.p,shape)
        q=np.broadcast_to(self.q,shape)
        self.temporal,self.spatial,self.S=matter_arrays(
            r,t,p,q,basis,shape)


def execute():
    rows=[]
    progress=HERE/'alpha_beta_c3_coupled_progress.json'
    for n in LEVELS:
        resolution=[n+1,n,n]
        progress.write_text(json.dumps(dict(status='IN_PROGRESS',rows=rows,
            active_resolution=resolution),indent=2)+'\n')
        print(json.dumps(dict(event='begin_fresh_coupled_solve',
            resolution=resolution,epsilon=EPS,omega=OMEGA)),flush=True)
        g=PhaseGrid(*resolution)
        u,w,it=g.solve()
        scalar=u.reshape(g.shape)
        vector=w.reshape((3,)+g.shape)
        audit=PhaseAudit(g.shape)
        regions=audit.summarize(scalar,vector)
        # Solver residuals use candidate-consistent S and E.
        diagnostics=g.diagnostics(u,w)
        inner=regions['inner_boundary_band']
        row=dict(resolution=resolution,coupled_iterations=int(it),
            solver_diagnostics=diagnostics,inner=inner,
            regions=regions)
        rows.append(row)
        np.savez_compressed(HERE/f'alpha_beta_c3_fields_{n+1}_{n}.npz',
            psi=scalar,W_cartesian=vector)
        print(json.dumps(dict(event='finished_fresh_solve',
            resolution=resolution,inner=inner)),flush=True)
    metrics=('independent_flux_H','independent_tensor_M_order2')
    orders={}
    for metric in metrics:
        values=[row['inner'][metric]['volume_weighted_rms'] for row in rows]
        orders[metric]=[float(math.log2(a/b)) for a,b in zip(values,values[1:])]
    gate={metric:all(p>=1.5 for p in values) for metric,values in orders.items()}
    result=dict(scope='EXPERIMENTAL changed matter phase source; NOT Run43 validation',
        parameters=dict(epsilon=EPS,omega=OMEGA,width=WIDTH,levels=LEVELS,
                        radial_switch='C3 septic'),
        fields='Theta=KTH*chi+epsilon*cos(chi+omega*t); Psi=KPS*chi+epsilon*s(r)*sin(chi-omega*t), at t=0',
        equations='Original conformal Hamiltonian and vector momentum solve; new spatial and temporal matter projections computed from the same scalar perturbations',
        rows=rows,orders=orders,inner_gate_experimental=gate,
        pass_experimental_both=all(gate.values()),
        note='New source differs from Run38/Run43; no claim of Run43 ninth check')
    (HERE/'alpha_beta_c3_coupled_results.json').write_text(json.dumps(result,indent=2)+'\n')
    progress.write_text(json.dumps(dict(status='COMPLETE',rows=rows),indent=2)+'\n')
    print(json.dumps(dict(status='COMPLETE',orders=orders,gate=gate)),flush=True)
    return int(not all(gate.values()))


if __name__=='__main__':
    raise SystemExit(execute())
