"""New radial energy certificate and resolution-dependent spectral bounds.

For even n>=8, let A=-4h^2 sym((D D)_II). Let S be the two independent
Dirichlet tridiagonal graph Laplacians on odd/even indices (edges at distance 2).
A-S/2 is a diagonally dominant symmetric M-matrix with nonnegative row sums,
hence PSD. Therefore -sym(L)>=S/(8h^2)>0. The longest chain has n/2 nodes,
giving lambda_min >= sin(pi/(n+2))^2/(2h^2), uniformly >=32 for Lshell=.2.
This proves stability of the RADIAL principal closure, not the full toroidal
variable-coefficient vector operator. Fourth-D2 spectra below are empirical.
"""
import sys,json
from pathlib import Path
import numpy as np
from run43_fourth_vector import stencil
from diagnose_toroidal_boundary_convergence import radial_matrix


def main():
    target=Path(sys.argv[1])
    if target.exists():raise FileExistsError(target)
    rows=[]
    for n in (8,16,32,64,128,256):
        h=.2/n;d=radial_matrix(n+1,h,2);L=(d@d)[1:-1,1:-1]
        A=-4*h*h*(L+L.T)/2
        S=2*np.eye(n-1)-np.eye(n-1,k=2)-np.eye(n-1,k=-2)
        remainder=A-S/2
        off=remainder-np.diag(np.diag(remainder))
        assert np.max(off)<1e-12 and np.min(np.sum(remainder,axis=1))>-1e-12
        lower=float(np.sin(np.pi/(n+2))**2/(2*h*h))
        observed=float(np.min(np.linalg.eigvalsh(-(L+L.T)/2)))
        assert observed>=lower-1e-8 and lower>=32-1e-8
        fourth=stencil(n+1,h,2).toarray()[1:-1,1:-1]
        rows.append(dict(n=n,radial_energy_lower_bound=lower,observed_energy_lower_bound=observed,original_singular_min=float(np.linalg.svd(L,compute_uv=False)[-1]),fourth_energy_lower_bound=float(np.min(np.linalg.eigvalsh(-(fourth+fourth.T)/2))),fourth_largest_eigenvalue_real=float(np.max(np.linalg.eigvals(fourth).real))))
    target.write_text(json.dumps(dict(scope='Radial principal-part stability only; not a full toroidal stability theorem.',uniform_lower_bound=32.,rows=rows),indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
