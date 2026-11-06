"""Conditional independent audit of four saved accepted grids; zero solves."""
from pathlib import Path
import sys, hashlib, json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'run43_final_validation'))
from evaluate import Evaluator, BASE, ROOT
OUT=Path(__file__).resolve().parent
def main():
    assert json.loads((OUT/'qualification.json').read_text())['status']=='PASS'
    assert json.loads((OUT/'continuum_identity_checks.json').read_text())['status']=='PASS'
    paths=[ROOT/'research/RUN43_CONSOLIDATED_MATHEMATICAL_SPECIFICATION_20261106.md',ROOT/'docs/RECOVERED_NONLINEAR_DEE_FORMULATION.md',ROOT/'reproduce/dee_stress_energy_3p1.py',ROOT/'reproduce/validate_toroidal_helical_dee_source.py',ROOT/'reproduce/solve_toroidal_helical_coupled_constraints.py',ROOT/'reproduce/diagnose_toroidal_boundary_convergence.py',ROOT/'research/audit_run43_arrays.py',ROOT/'research/RUN43_COORDINATED_REPLAY_20261103.md']
    paths+=list(BASE.glob('*')); paths+=[ROOT/'research/results/coordinated-replay-20261103-1831/coordination/frozen_audit.py',ROOT/'research/results/coordinated-replay-20261103-1831/validation/baseline-independent-validation.json']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
    (OUT/'input_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    result={'qualification_status':'PASS: see qualification.json; rates resolved above floating-point bound','components':['H','Mx','My','Mz'],'reference_commit':'5799e440b5f27bfe81c3ce5c44814a7d2d652157','rows':[],'orders':[]}
    for n in (8,16,32,64):
        path=BASE/f'fields_{n}.npz';expected=json.loads((BASE/f'result_{n}.json').read_text())['sha256'];assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
        with np.load(path) as data:u=data['psi'];w=data['W_cartesian']
        assert u.shape==(n+1,n,n) and w.shape==(3,n+1,n,n)
        assert np.isfinite(u).all() and np.isfinite(w).all() and u.min()>0
        e=Evaluator(n,9); high,scale=e.residual(u,w)
        mid,_=Evaluator(n,7).residual(u,w);low,_=Evaluator(n,5).residual(u,w)
        angular,_=Evaluator(n,9,n//3).residual(u,w)
        uncertainty=np.maximum.reduce([abs(high-mid),abs(mid-low),abs(high-angular)])
        boundary={'psi_trace_max_error':float(np.max(abs(u[[0,-1]]-1))),'W_trace_max_error':float(np.max(abs(w[:,[0,-1]]))),'angular_seams':'Endpoint-exclusive arrays: no independent duplicate 2pi trace is saved; equality of physical seam traces NOT YET ESTABLISHED. Fourier continuation is periodic by construction, not evidence of saved endpoint equality.','adjacent_wrap_max':{}}
        for axis in (1,2):
            boundary['adjacent_wrap_max'][str(axis)]={'psi':float(np.max(abs(np.take(u,0,axis)-np.take(u,-1,axis)))),'W':float(np.max(abs(np.take(w,0,axis+1)-np.take(w,-1,axis+1))))}
        row={'n':n,'boundary':boundary,'regions':{}}
        for name,ids in e.regions().items():
            v=e.rms(high,u,ids);err=e.rms(uncertainty,u,ids);den=e.rms(scale,u,ids)
            row['regions'][name]={'rms':v.tolist(),'max_abs':np.max(abs(high[:,ids]),axis=(1,2,3)).tolist(),'term_scale_rms':den.tolist(),'normalized_rms':(v/den).tolist(),'uncertainty_rms':err.tolist(),'uncertainty_fraction':(err/v).tolist(),'within_10_percent':(err<=.1*v).tolist()}
        result['rows'].append(row)
        print('evaluated',n,flush=True)
    for name in result['rows'][0]['regions']:
        for a,b in zip(result['rows'],result['rows'][1:]):
            av=a['regions'][name];bv=b['regions'][name]
            aa=np.array(av['rms']);bb=np.array(bv['rms']);ua=np.array(av['uncertainty_rms']);ub=np.array(bv['uncertainty_rms'])
            p=np.log2(aa/bb)
            bound=np.full(4,-np.inf); valid=aa>ua
            bound[valid]=np.log2((aa[valid]-ua[valid])/(bb[valid]+ub[valid]))
            result['orders'].append({'region':name,'pair':[a['n'],b['n']],'nominal':p.tolist(),'conservative':[float(v) if np.isfinite(v) else None for v in bound],'passes_p_1_5':(bound>=1.5).tolist(),'both_uncertainties_within_budget':(np.maximum(ua/aa,ub/bb)<=.1).tolist()})
    result['uncertainty_scope']='Empirical envelope of 5/7/9-point radial reconstructions and removal of angular modes above n//3. Null conservative order means uncertainty contains zero; no finite lower order established. Not a rigorous error bound; unresolved angular aliases cannot be bounded from saved samples alone.'
    result['absolute_accuracy']='NOT YET ESTABLISHED: no predeclared absolute physical residual tolerance recovered; algebraic solver tolerances are not physical tolerances.'
    (OUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
