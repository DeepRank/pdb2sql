from time import time
import numpy as np
from pdb2sql.pdb2sqlAlchemy import pdb2sql_alchemy
import unittest


class TestAlchemy(unittest.TestCase):
    """Test Alchemy generation."""

    def test_alchemy(self):

        # create the data base
        t0 = time()
        db = pdb2sql_alchemy(self.pdb)
        print('ALCH %f' % (time() - t0))

        # extract the xyz position
        xyz = db.get('x,y,z', model=0)

        # put the data back
        db.update('x,y,z', xyz)

        # extract atoms but not 'H' atoms
        atoms = db.get(chainID='A',resName=['THR','GLN'],no_element=['H'], model=0)
        for at in atoms:
            print(at.resName, at.name, at.x, at.y, at.z, at.element)

    def setUp(self):
        self.pdb = 'pdb/test_model.pdb'


if __name__ == '__main__':
    unittest.main()
