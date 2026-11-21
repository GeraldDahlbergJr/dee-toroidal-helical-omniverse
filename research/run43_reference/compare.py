"""Use the frozen audit operator for a weighted vector residual decomposition."""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
from pathlib import Path

import numpy as np
import scipy
from reference import ROOT
from diagnose_toroidal_boundary_convergence import Audit


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_cartesian(path, n, nr):
    with np.load(path) as data:
        u = data['psi'][::nr//n]
        wr,wp,wz = data['W'][:,::nr//n,: ,None]
    phi = np.arange(n)*2*np.pi/n
    w = np.array(np.broadcast_arrays(wr*np.cos(phi)-wp*np.sin(phi),wr*np.sin(phi)+wp*np.cos(phi),wz))
    return np.broadcast_to(u[:,:,None],(n+1,n,n)),w


def csv_table(path, rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def analyze(data,hybrid,out):
    out.mkdir(parents=True)
    native = out/'native'
    native.mkdir()
    previous = json.loads((ROOT/'research/run43_hybrid/measurements.json').read_text())
    for path,sha in previous['sha256'].items():
        assert digest(ROOT/path)==sha, path
    bands, rows, refinements, histories = [],[],[],[]
    for i,n in enumerate((8,16,32,64)):
        hpath = hybrid/f'fields_{n}.npz'
        old = previous['experiments']['hybrid']['levels'][i]
        assert digest(hpath)==old['sha256']
        with np.load(hpath) as z:
            u,w = z['psi'],z['W_cartesian']
        audit = Audit(u.shape)
        rv = audit.r[:,0,0]
        indices = np.flatnonzero((rv<=.10+1e-12)&(np.arange(n+1)>0))
        mass = (audit.r*audit.q*u**6)[indices]
        mh = audit.momentum(u,w,2)[:,indices]
        def inner(x,y):
            return (mass*(x*y).sum(axis=0)).sum(axis=(1,2))
        prev = None
        for factor in (1,2,4,8):
            nr=n*factor
            folder=data/f'n{n}-r{nr}'
            hist=json.loads((folder/'result.json').read_text())
            assert digest(folder/'reduced.npz')==hist['sha256']
            histories.append(hist)
            ur,wr=load_cartesian(folder/'reduced.npz',n,nr)
            mr=audit.momentum(ur,wr,2)[:,indices]
            if factor==1:
                path=native/f'fields_{n}.npz'
                np.savez_compressed(path,psi=ur,W_cartesian=wr)
                (native/f'result_{n}.json').write_text(json.dumps(dict(sha256=digest(path))))
            diff = mh-mr
            numbers = dict(reference=inner(mr,mr),difference=inner(diff,diff),
                           cross=2*inner(mr,diff),hybrid=inner(mh,mh))
            massrow=mass.sum(axis=(1,2))
            def summarize(index):
                m=massrow[index].sum()
                vals={k:float(v[index].sum()/m) for k,v in numbers.items()}
                closure=vals['hybrid']-vals['reference']-vals['difference']-vals['cross']
                assert abs(closure)<1e-12*vals['hybrid']
                vals.update(reference_rms=np.sqrt(vals['reference']),difference_rms=np.sqrt(vals['difference']),
                            hybrid_rms=np.sqrt(vals['hybrid']),closure=closure)
                if prev is not None:
                    change=mr-prev[0]
                    vals['last_ref_change_rms']=float(np.sqrt(inner(change,change)[index].sum()/m))
                return vals
            band = summarize(slice(None))
            assert abs(band['hybrid_rms']-old['band_rms']['inner_boundary_band']['independent_tensor_M_order2'])<1e-15
            refinements.append(dict(n=n,nr=nr,ref_rms=band['reference_rms'],
                ref_vector_change=band.get('last_ref_change_rms'),max_psi_change=None if prev is None else float(abs(ur-prev[1]).max()),
                max_W_change=None if prev is None else float(abs(wr-prev[2]).max()),
                mms_max=hist['mms']['max_node_error'],mms_first=hist['mms']['first_audit_row_error']))
            if factor==8:
                bands.append(dict(n=n,nr=nr,**band))
                for k,j in enumerate(indices):
                    rowvals=summarize(k)
                    fraction=massrow[k]/massrow.sum()
                    rows.append(dict(n=n,row=int(j),radius=float(rv[j]),band_mass_fraction=fraction,
                                     **rowvals,**{name+'_contribution_to_band':fraction*rowvals[name]
                                     for name in ('reference','difference','cross','hybrid')}))
            prev = mr,ur,wr
    csv_table(out/'bands.csv',bands)
    csv_table(out/'rows.csv',rows)
    csv_table(out/'refinement.csv',refinements)
    pairs=[]
    for a,b in zip(bands,bands[1:]):
        def energy_change(level):
            eps=level['last_ref_change_rms']
            return 2*level['reference_rms']*eps+eps**2
        pairs.append(dict(last_ref_deficit_sensitivity=energy_change(b)+energy_change(a)/8,
            coarse=a['n'],fine=b['n'],
            p_hybrid=np.log2(a['hybrid_rms']/b['hybrid_rms']),
            p_ref=np.log2(a['reference_rms']/b['reference_rms']),
            p_delta=np.log2(a['difference_rms']/b['difference_rms']),
            **{k+'_gate_deficit':b[k]-a[k]/8 for k in ['reference','difference','cross','hybrid']}))
    csv_table(out/'pairs.csv',pairs)
    (out/'histories.json').write_text(json.dumps(histories,indent=2)+'\n')
    paths=list((ROOT/'research/run43_reference').glob('*.py'))+[ROOT/p for p in previous['sha256']]
    (out/'provenance.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,
        scipy=scipy.__version__,sha256={str(p.relative_to(ROOT)):digest(p) for p in paths},
        hybrid_sha256=[v['sha256'] for v in previous['experiments']['hybrid']['levels']]),indent=2)+'\n')
    audit=subprocess.run([__import__('sys').executable,str(ROOT/'research/audit_run43_arrays.py'),str(native)],
                         stdout=subprocess.PIPE,text=True)
    (out/'native-audit.log').write_text(audit.stdout)
    report=json.loads((native/'independent-audit.json').read_text())
    assert len(report['checks'])==12 and audit.returncode==(0 if report['passed']==12 else 1)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    for name in ['data','hybrid','out']:
        cli.add_argument(name,type=Path)
    a=cli.parse_args()
    analyze(a.data,a.hybrid,a.out)
