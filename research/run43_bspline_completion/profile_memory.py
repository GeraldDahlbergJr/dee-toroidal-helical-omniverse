"""Record n=32 forcing and individual weak applications, without solving."""
import json
import os
import resource
import time

os.sched_setaffinity(0,{max(os.sched_getaffinity(0))})
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))

import experiment as ex
import numpy as np

started=time.monotonic()
profile=[]
def checkpoint(phase):
    profile.append({
        'phase':phase,
        'elapsed_s':time.monotonic()-started,
        'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    })
    print(json.dumps(profile[-1]),flush=True)

space=ex.Space(32,6)
checkpoint('space')
forcing_blocks=[ex.forcing(block) for block in space.slabs()]
checkpoint('forcing')
rng=np.random.default_rng(4332)
c=rng.normal(scale=.001,size=space.cshape)
w=rng.normal(scale=.0001,size=(3,)+space.cshape)
space.vector(w)
checkpoint('one_vector_application')
ex.vector_rhs(space,c,forcing_blocks)
checkpoint('one_vector_load')
sq_blocks=[ex.strain_squared(block,w) for block in space.slabs()]
checkpoint('strain_squared')
_,potential=ex.scalar_rhs(space,c,forcing_blocks,sq_blocks)
checkpoint('one_scalar_load')
space.jacobian_scalar(c,potential)
checkpoint('one_scalar_application')
