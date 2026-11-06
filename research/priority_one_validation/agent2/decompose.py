"""Read-only decomposition of archived evaluator sensitivity; no solves."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/run43_final_validation'))
from evaluate import Evaluator,BASE
out={'purpose':'Decompose existing envelope; not a corrected evaluator or physical acceptance','rows':[]}
for n in (8,16,32,64):
 p=BASE/f'fields_{n}.npz'
 assert hashlib.sha256(p.read_bytes()).hexdigest()==json.loads((BASE/f'result_{n}.json').read_text())['sha256']
 with np.load(p) as z:u=z['psi'];w=z['W_cartesian']
 e=Evaluator(n,9);h,_=e.residual(u,w)
 mid,_=Evaluator(n,7).residual(u,w);low,_=Evaluator(n,5).residual(u,w)
 cut,_=Evaluator(n,9,n//3).residual(u,w)
 diffs={'radial_9_7':abs(h-mid),'radial_7_5':abs(mid-low),'angular_cut':abs(h-cut)}
 envelope=np.maximum.reduce(list(diffs.values()))
 row={'n':n,'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'regions':{},'spectral':{}}
 for region,ids in e.regions().items():
  r=e.rms(h,u,ids)
  row['regions'][region]={'residual_rms':r.tolist(),'envelope_fraction':(e.rms(envelope,u,ids)/r).tolist(),'difference_fraction':{k:(e.rms(v,u,ids)/r).tolist() for k,v in diffs.items()}}
 for name,field in [('psi_minus_1',u-1),('W',w)]:
  ax1=field.ndim-2;ax2=field.ndim-1
  f=np.fft.fftn(field,axes=(ax1,ax2));power=abs(f)**2
  k=np.fft.fftfreq(n,1/n);mask=(abs(k[:,None])>n//3)|(abs(k[None,:])>n//3)
  row['spectral'][name]={'fraction_power_above_cut':float(power[...,mask].sum()/power.sum()),'fraction_power_theta_nyquist':float(power[...,n//2,:].sum()/power.sum()),'fraction_power_phi_nyquist':float(power[...,:,n//2].sum()/power.sum())}
 wt=np.broadcast_to(e.r*e.q,e.shape)*u**6
 energy=np.sum((h-mid)**2*wt,axis=(2,3))
 total=np.sum(energy[:,1:-1],axis=1)
 row['radial_9_7_sensitivity_energy']={'per_interior_row_fraction':(energy[:,1:-1]/total[:,None]).tolist(),'first_four_fraction':(energy[:,1:5].sum(axis=1)/total).tolist(),'last_four_fraction':(energy[:,-5:-1].sum(axis=1)/total).tolist()}
 prior=json.loads((ROOT/'research/run43_nonpolynomial_validation/results.json').read_text())['rows']
 row['existing_envelope_max_relative_difference']=max(float(np.max(abs(np.array(row['regions'][name]['envelope_fraction'])/np.array(prior[len(out['rows'])]['regions'][name]['uncertainty_fraction'])-1))) for name in row['regions'])
 out['rows'].append(row)
 print(n,flush=True)
(Path(__file__).parent/'decomposition.json').write_text(json.dumps(out,indent=2)+'\n')
