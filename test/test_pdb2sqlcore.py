from time import time
import numpy as np
from pdb2sql import pdb2sql
import unittest


class TestCore(unittest.TestCase):
    """Test Core generation."""

    def test_core(self):

        t0 = time()
        db = pdb2sql(self.pdb)
        print('SQL %f' % (time() - t0))
        xyz = db.get('x,y,z', model=0)
        # db.update('x,y,z',xyz)

        xyz = db.get('x,y,z',
            chainID='A',
            resName=['VAL', 'LEU'],
            name=[ 'CA', 'C', 'O', 'N'])

        xyz = db.get('x,y,z', chainID='A', resSeq=1)
        xyz = np.array(xyz)
        xyz -= np.mean(xyz)
        # db.update('x,y,z',xyz,chainID='A',resSeq=1)

    def setUp(self):
        self.pdb = 'pdb/5hvd.pdb'


if __name__ == '__main__':
    unittest.main()
