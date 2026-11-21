"""Check radial discrete Green identity and physical source at quadrature points."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse as sp
from reference import Reference, quadrature, RHO, potential
from validate_toroidal_helical_dee_source import source


def run():
    results=[]
    for n in (8,16,32,64):
        ref=Reference(n,n)
        p,d,r,w=quadrature(n//4)
        trace=np.zeros(n+1)
        trace[0],trace[-1]=-1,1
        defect=p.T@sp.diags(w)@d+d.T@sp.diags(w)@p-sp.diags(trace)
        errors=[]
        for k in np.linspace(0,len(ref.r)-1,40,dtype=int):
            theta=2*np.pi*(k%n)/n
            r=ref.r[k]
            energy,cov,_,_=source(r,theta)
            s,c=np.sin(theta),np.cos(theta)
            cart=np.array([cov[0]*c-cov[1]*s/r,cov[2]/ref.q[k],cov[0]*s+cov[1]*c/r])
            errors.append(max(abs(cart-ref.source[:,k]).max(),
                          abs(energy-(.5*(ref.temporal+ref.spatial[k])+potential(RHO)))))
        results.append(dict(n=n,green_max=float(abs(defect).max()),source_max=float(max(errors)),
                            stiffness_skew_max=float(abs(ref.K-ref.K.T).max())))
    assert max(max(row['green_max'],row['source_max'],row['stiffness_skew_max']) for row in results)<1e-11
    return results


if __name__=='__main__':
    Path(sys.argv[1]).write_text(json.dumps(run(),indent=2)+'\n')
