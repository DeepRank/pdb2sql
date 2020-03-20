import os 
import unittest
from pdb2sql.align import align, pca
import numpy as np 

class TestAlign(unittest.TestCase):
    """Test the superpose functionality"""

    def setUp(self):
        self.pdb = 'pdb/1AK4/1AK4_10w.pdb' 

    def test_align(self):
        """test the routine."""

        for idir, axis in zip([0,1,2],['x','y','z']):
            sql = align(self.pdb, axis=axis)
            xyz = np.array(sql.get('x,y,z'))
            u,v,_ = pca(xyz)
            vmax = np.abs(v[:,np.argmax(u)])
            assert np.argmax(vmax) == idir


