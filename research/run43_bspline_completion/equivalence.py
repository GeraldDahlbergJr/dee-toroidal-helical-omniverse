"""Compare streamed and historical forcing and weak kernels without PDE solves."""
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time

os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
import experiment as fixed
import numpy as np

historical_file=Path(sys.argv[1])
spec=importlib.util.spec_from_file_location('historical_experiment',historical_file)
historical=importlib.util.module_from_spec(spec)
spec.loader.exec_module(historical)
del fixed.FUN
fixed.FUN=historical.FUN

rtol=5e-12
atol=5e-14
rows=[]
for n in (8,16):
    start=time.monotonic()
    old=historical.Space(n,6)
    new=fixed.Space(n,6)
    rng=np.random.default_rng(4316)
    c=rng.normal(scale=.001,size=old.cshape)
    w=rng.normal(scale=.0001,size=(3,)+old.cshape)
    old_forcing=historical.forcing(old)
    new_forcing=[fixed.forcing(block) for block in new.slabs()]
    results={}
    def check(name,a,b):
        np.testing.assert_allclose(a,b,rtol=rtol,atol=atol)
        results[name]={
            'max_abs':float(np.max(np.abs(a-b))),
            'relative_l2':float(np.linalg.norm((a-b).ravel())/np.linalg.norm(a.ravel())),
        }

    check('forcing_scalar',old_forcing[0],np.concatenate([v[0] for v in new_forcing],axis=0))
    check('forcing_vector',old_forcing[1],np.concatenate([v[1] for v in new_forcing],axis=1))
    check('scalar_principal',old.scalar(c),new.scalar(c))
    check('vector_principal',old.vector(w),new.vector(w))
    u=1+old.val(c)
    sc,vc=old.sources(u)
    check('vector_load',np.array([old.back(old_forcing[1][i]-vc[i]) for i in range(3)]),
          fixed.vector_rhs(new,c,new_forcing))
    dw=np.array([old.grad(z) for z in w])
    div=sum(dw[i,i] for i in range(3))
    strain=dw+dw.swapaxes(0,1)
    for i in range(3):
        strain[i,i]-=2*div/3
    sq=np.sum(strain**2,axis=(0,1))
    sq_blocks=[fixed.strain_squared(block,w) for block in new.slabs()]
    old_rhs=old.back(old_forcing[0]+sq*u**-7/8+sc)-old.scalar(c)
    new_rhs,new_pot=fixed.scalar_rhs(new,c,new_forcing,sq_blocks)
    check('scalar_residual',old_rhs,new_rhs)
    old_bb=.068614375*(old.r**-2+9*old.q**-2)
    old_pot=-7*sq*u**-8/8+2*np.pi*historical.G*(5*.0127638*u**4+old_bb)
    check('scalar_jacobian',old.scalar(c)-old.back(old_pot*old.val(c)),new.jacobian_scalar(c,new_pot))
    rows.append({'n':n,'metrics':results,'seconds':time.monotonic()-start,
                 'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
print(json.dumps({'atol':atol,'rtol':rtol,'rows':rows},indent=2))
