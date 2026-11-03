import json,os,sys,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from run43_compatible_closure import CompatibleGrid
from run43_second_hessian import SecondHessianGrid
from check_run43_compatible_accuracy import check


class CompatibleBoundary(unittest.TestCase):
    def test_original_sources_and_analytic_second_order_consistency(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'flux','DEE_PRECONDITIONER':'fourier'}):
            for cls in (CompatibleGrid,SecondHessianGrid):
                report=check(cls)
                for region,orders in report['orders'].items():
                    self.assertGreater(orders[-1],1.5,(cls.__name__,region,orders))

    def test_small_grid_spectrum_has_no_zero_or_growing_eigenmode(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'2','DEE_PRECONDITIONER':'ilu'}):
            for cls in (CompatibleGrid,SecondHessianGrid):
                g=cls(6);identity=np.eye(g.vector.shape[0])
                matrix=np.column_stack([g.vector@x for x in identity])
                eig=np.linalg.eigvals(matrix)
                self.assertLess(np.max(eig.real),-1.,(cls.__name__,np.max(eig.real)))


if __name__=='__main__':unittest.main()
