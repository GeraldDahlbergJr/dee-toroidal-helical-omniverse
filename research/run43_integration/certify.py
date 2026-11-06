"""Enclose the exact tensor mass on ONE curved element; no operator assembly."""
import argparse
import json
import platform
from pathlib import Path
import flint
from flint import arb, arb_mat, ctx
from itertools import product
ALPHAS=list(product(range(3),repeat=3))
SELECTED=list(range(3,101))

def tensors(x):
 r=arb('1/10')+x[0]/100;t=x[1]/100;p=x[2]/100
 ct,st,cp,sp=t.cos(),t.sin(),p.cos(),p.sin();q=1+r*ct
 er=[ct*cp,ct*sp,st];et=[-st*cp,-st*sp,ct];ep=[-sp,cp,arb(0)]
 pairs=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
 ts=[]
 for a in ALPHAS:
  ds=[arb(a[k])*x[k]**max(0,a[k]-1)*x[(k+1)%3]**a[(k+1)%3]*x[(k+2)%3]**a[(k+2)%3] if a[k] else arb(0) for k in range(3)]
  g=[100*(er[j]*ds[0]+et[j]*ds[1]/r+ep[j]*ds[2]/q) for j in range(3)]
  for i in range(3):
   ts.append([(g[k] if j==i else 0)+(g[j] if k==i else 0)-(arb(2)*g[i]/3 if j==k else 0) for j,k in pairs])
 E=[[1,-1,0,0,0,0],[1,1,-2,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]
 for e in E:
  for v in [arb(1),q*cp,q*sp,r*st]:ts.append([v*c for c in e])
 return ts[3:],r*q/1000000

def sample_certificate():
 # 64 interior rational sample points; full Frobenius norm, no integration.
 rows=[]
 for ijk in product(range(1,5),repeat=3):
  ts,_=tensors([arb(i)/5 for i in ijk])
  for c in range(6):rows.append([v[c]*(arb(2).sqrt() if c>=3 else 1) for v in ts])
 S=arb_mat(rows);A=S.transpose()*S
 inverse=A.inv()
 # total absolute entry sum bounds spectral norm and avoids nonrigorous max.
 bound=sum((abs(inverse[i,j]) for i in range(98) for j in range(98)),arb(0))
 k=0
 while not bound < arb(2)**k:
  k+=1
  if k>4096:raise ArithmeticError('bounded inverse-norm search failed')
 a=(k+7)//2;d=a+33;e=23+3*d+2*a
 assert arb(2)**(-d) < arb(1)/5
 assert 3*arb(2)**30*arb(2)**(-d) < arb(2)**(-a)/2
 return dict(sample_count=64,sample_inverse_norm_upper_power2=k,alpha_lower_power2=-a,delta_power2=-d,mass_eigenvalue_lower_power2=-e)

def intervals(n):
 H=[[arb(0) for j in range(98)] for i in range(98)]
 for indices in product(range(n),repeat=3):
  xs=[arb(2*i+1)/(2*n)+arb(0,(arb(1)/(2*n)).upper()) for i in indices]
  ts,w=tensors(xs);w=w/n**3
  for i in range(98):
   for j in range(i,98):
    v=sum((ts[i][c]*ts[j][c]*(2 if c>=3 else 1) for c in range(6)),arb(0))
    H[i][j]+=w*v
 for i in range(98):
  for j in range(i):H[i][j]=H[j][i]
 return H
def dyadic(x):
 return list(map(int,x.man_exp()))

def test(H):
 # Independent closed-form volume check: E1:E1=2, E2:E2=6.
 r0,r1=arb('1/10'),arb('11/100')
 volume=arb('1/100')*(arb('1/100')*(r1*r1-r0*r0)/2
             +(r1**3-r0**3)*(arb('1/100').sin())/3)
 assert H[78][78].contains(2*volume)
 assert H[82][82].contains(6*volume)
 assert H[78][82].contains(0)
 assert H[86][86].contains(2*volume) # xy pair counts twice.
 assert all(dyadic(H[i][j].mid())==dyadic(H[j][i].mid())
            and dyadic(H[i][j].rad())==dyadic(H[j][i].rad())
            for i in range(98) for j in range(98))
 ts,w=tensors([arb(1)/2]*3)
 assert w > arb(2)**-24
 assert len(ts)==98 and all((v[0]+v[1]+v[2]).contains(0) for v in ts)
 return ['closed_form_volume_E1_E2_E3','orthogonal_STF_pair',
         'identical_symmetric_enclosures','positive_Jacobian','trace_free']

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--subdivisions',type=int,default=4)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 if not 1<=args.subdivisions<=8:parser.error('bounded test requires 1..8 subdivisions')
 ctx.prec=1024
 certificate=sample_certificate()
 ctx.prec=128
 H=intervals(args.subdivisions)
 tests=test(H)
 entries=[]
 for i in range(98):
  for j in range(i,98):
   b=H[i][j]
   assert b.is_finite()
   entries.append([i,j,dyadic(b.mid()),dyadic(b.rad())])
 error_power=0
 while not all(b.rad()<arb(2)**error_power for row in H for b in row):error_power+=1
 while all(b.rad()<arb(2)**(error_power-1) for row in H for b in row):error_power-=1
 result=dict(checkpoint='ade8206d0f629fd8a73e6f66ef6490f0b1f5084d',
   method='validated interval subdivision (exact target, enclosure output)',
   python=platform.python_version(),python_flint=flint.__version__,
   flint=flint.__FLINT_VERSION__,sample_precision_bits=1024,
   integration_precision_bits=128,subdivisions_per_axis=args.subdivisions,
   reference_element={'r':['1/10','11/100'],'theta':['0','1/100'],'phi':['0','1/100']},
   raw_generator_count=101,discarded_zero_generators=[0,1,2],
   selected_generators=SELECTED,dimension=98,certificate=certificate,
   all_entry_absolute_errors_strictly_below_power2=error_power,
   tests=tests,upper_triangle_entry_format=['i','j','midpoint_[integer,exponent]',
                                            'radius_[integer,exponent]'],
   upper_triangle=entries)
 args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items()
                   if k not in ('upper_triangle','selected_generators')},indent=2))

if __name__=='__main__':main()
