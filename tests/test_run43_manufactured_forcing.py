"""Check the new analytic forcing by independent Cartesian differentiation."""
import sys,os,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from run43_manufactured_boundary_solve import analytic
from run43_fourth_vector import FourthVectorGrid


class AnalyticForcing(unittest.TestCase):
    def test_forcing_against_50_digit_cartesian_derivatives(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'flux','DEE_PRECONDITIONER':'fourier'}):g=FourthVectorGrid(8)
        with mp.workdps(50):
            c=[mp.mpf('.3'),mp.mpf('-.4'),mp.mpf('.5')]
            for kind in ('radial','mixed'):
                _,truth=analytic(g,kind)
                for ijk in ((1,1,1),(3,2,2),(5,3,5)):
                    k=np.ravel_multi_index(ijk,g.shape)
                    x,y,z=[mp.mpf(str(v)) for v in (g.q[k]*np.cos(g.p[k]),g.q[k]*np.sin(g.p[k]),g.r[k]*np.sin(g.t[k]))]
                    def F(x,y,z):
                        radius=mp.sqrt((mp.sqrt(x*x+y*y)-1)**2+z*z)
                        f=mp.mpf('.30')-radius-mp.mpf('.0125')/radius
                        return f if kind=='radial' else f*z/radius
                    hess=[]
                    for i in range(3):
                        row=[]
                        for j in range(3):
                            degree=[0,0,0];degree[i]+=1;degree[j]+=1
                            row.append(mp.diff(F,(x,y,z),tuple(degree)))
                        hess.append(row)
                    lap=sum(hess[i][i] for i in range(3))
                    expected=[float(c[i]*lap+sum(hess[i][j]*c[j] for j in range(3))/3) for i in range(3)]
                    np.testing.assert_allclose(truth[:,k],expected,rtol=2e-12,atol=2e-11)


if __name__=='__main__':unittest.main()
