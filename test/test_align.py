import unittest
from pathlib import Path
from pdb2sql.align import align, pca, align_interface
import numpy as np


from . import pdb_folder


class TestAlign(unittest.TestCase):
    """Test the superpose functionality"""

    def setUp(self):
        self.pdb = Path(pdb_folder, '1AK4', '1AK4_10w.pdb')

    def test_align(self):
        """Test align()"""

        for idir, axis in zip([0, 1, 2], ['x', 'y', 'z']):
            with self.assertWarns(UserWarning) as ex:
                db = align(self.pdb, axis=axis)
            xyz = np.array(db.get('x,y,z'))
            u, v = pca(xyz)
            vmax = np.abs(v[:, np.argmax(u)])
            assert np.argmax(vmax) == idir

    def test_align_interface(self):
        """Test align_interface()"""

        def get_xyz_interface(db):
            idx = db.get_contact_atoms()
            idx = idx['A'] + idx['B']
            return np.array(db.get('x,y,z', rowID=idx))

        for idir, plane in zip([2, 1, 0], ['xy', 'xz', 'yz']):
            with self.assertWarns(UserWarning) as ex:
                db = align_interface(self.pdb, plane=plane)
            xyz = get_xyz_interface(db)
            u, v = pca(xyz)
            vmax = np.abs(v[:, np.argmin(u)])
            assert np.argmax(vmax) == idir

if __name__ == '__main__':
    unittest.main()
