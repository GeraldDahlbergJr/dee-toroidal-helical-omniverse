import os,sys,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from compare_run43_operators import radial_errors
from run43_endpoint_continuation import EndpointGrid
from run43_consistent_momentum import TensorGrid


class EndpointContinuation(unittest.TestCase):
    def test_matched_composition_first_interior_rows_converge(self):
        rows=radial_errors()
        for region in ('first_inner','first_outer','endpoint_max'):
            p=[np.log2(a['matched'][region]/b['matched'][region]) for a,b in zip(rows,rows[1:])]
            self.assertTrue(all(x>1.5 for x in p),p)
        # Full vector endpoint PDE compatibility remains a separate diagnostic.

    def test_matter_is_unchanged_under_endpoint_correction(self):
        with patch.dict(os.environ,{'DEE_SCALAR_ORDER':'flux','DEE_PRECONDITIONER':'fourier'}):
            g=TensorGrid(8);c=EndpointGrid(8)
            u=1+np.linspace(0,.001,g.size)
            np.testing.assert_array_equal(c.S,g.S)
            np.testing.assert_array_equal(c.energy(u),g.energy(u))
            np.testing.assert_array_equal(c.L.toarray(),g.L.toarray())

if __name__=='__main__':unittest.main()
