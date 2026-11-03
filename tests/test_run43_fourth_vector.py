import os,sys,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from run43_fourth_vector import stencil,FourthVectorGrid
from run43_consistent_momentum import TensorGrid
from check_run43_fourth_accuracy import check


class FourthVector(unittest.TestCase):
    def test_radial_stencil_polynomial_exactness_including_endpoints(self):
        x=np.linspace(.05,.25,9);h=x[1]-x[0]
        for degree in (1,2):
            d=stencil(len(x),h,degree)
            for power in range(5 if degree==1 else 6):
                exact=np.zeros_like(x) if power<degree else (power if degree==1 else power*(power-1))*x**(power-degree)
                np.testing.assert_allclose(d@x**power,exact,atol=5e-11,rtol=1e-9)

    def test_analytic_radial_and_angular_accuracy(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'flux','DEE_PRECONDITIONER':'fourier'}):r=check((8,16,32))
        for region,orders in r['radial_orders'].items():self.assertTrue(all(p>3 for p in orders),(region,orders))
        self.assertTrue(all(p>3 for p in r['cartesian_polynomial_orders']),r)

    def test_original_matter_and_scalar_laplacian_preserved(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'flux','DEE_PRECONDITIONER':'fourier'}):g=FourthVectorGrid(8);b=TensorGrid(8)
        np.testing.assert_array_equal(g.S,b.S)
        u=np.linspace(1,1.01,g.size);np.testing.assert_array_equal(g.energy(u),b.energy(u))
        np.testing.assert_array_equal(g.lap.toarray(),b.lap.toarray())

    def test_small_grid_has_no_exact_null_mode(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'2','DEE_PRECONDITIONER':'ilu'}):g=FourthVectorGrid(6)
        matrix=np.column_stack([g.vector@x for x in np.eye(g.vector.shape[0])]);s=np.linalg.svd(matrix,compute_uv=False)
        self.assertGreater(s[-1]/s[0],1e-7)

if __name__=='__main__':unittest.main()
