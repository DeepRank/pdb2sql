import unittest
from pdb2sql.align import align, pca, align_interface
import numpy as np


class TestAlign(unittest.TestCase):
    """Test the superpose functionality"""

    def setUp(self):
        self.pdb = 'pdb/1AK4/1AK4_10w.pdb'

    def test_align(self):
        """test the routine."""

        for idir, axis in zip([0, 1, 2], ['x', 'y', 'z']):
            sql = align(self.pdb, axis=axis)
            xyz = np.array(sql.get('x,y,z'))
            u, v = pca(xyz)
            vmax = np.abs(v[:, np.argmax(u)])
            assert np.argmax(vmax) == idir

    def test_align_interface(self):
        """test the routine."""

        def get_xyz_interface(sql):
            idx = sql.get_contact_atoms()
            idx = idx['A'] + idx['B']
            return np.array(sql.get('x,y,z', rowID=idx)) 

        for idir, plane in zip([2, 1, 0], ['xy', 'xz', 'yz']):
            sql = align_interface(self.pdb, plane=plane)
            xyz = get_xyz_interface(sql)
            u, v = pca(xyz)
            vmax = np.abs(v[:, np.argmin(u)])
            assert np.argmax(vmax) == idir

    


