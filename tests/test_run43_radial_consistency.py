import sys,os,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from run43_radial_consistency import RadialGrid
from run43_consistent_momentum import TensorGrid
from check_run43_radial_accuracy import check


class RadialConsistency(unittest.TestCase):
    def test_radial_principal_symbol_correction_has_no_first_order_edge_defect(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'2','DEE_PRECONDITIONER':'ilu'}):
            report=check((8,16,32))
        for region,orders in report['orders']['corrected'].items():
            self.assertTrue(all(p>1.5 for p in orders),(region,orders))

    def test_original_sources_and_scalar_operator_preserved(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'flux','DEE_PRECONDITIONER':'fourier'}):
            g=RadialGrid(8);b=TensorGrid(8)
        np.testing.assert_array_equal(g.S,b.S)
        np.testing.assert_array_equal(g.energy(np.linspace(1,1.01,g.size)),b.energy(np.linspace(1,1.01,b.size)))
        np.testing.assert_array_equal(g.lap.toarray(),b.lap.toarray())

    def test_radial_checkerboard_damping_and_small_grid_invertibility(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'2','DEE_PRECONDITIONER':'ilu'}):g=RadialGrid(6)
        eye=np.eye(g.vector.shape[0]);matrix=np.column_stack([g.vector@x for x in eye])
        s=np.linalg.svd(matrix,compute_uv=False)
        self.assertGreater(s[-1]/s[0],1e-7)
        pattern=(-1.)**np.indices(g.shape)[0].ravel();pattern.reshape(g.shape)[[0,-1]]=0
        w=np.array([pattern,np.zeros(g.size),np.zeros(g.size)])
        x=w[:,g.idx].ravel()
        self.assertLess(float(x@(g.vector@x)),0.)

if __name__=='__main__':unittest.main()
