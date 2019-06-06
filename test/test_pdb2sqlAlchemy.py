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
        print('ALCH %f' %(time()-t0))

        # extract the xyz position of all VAL and LEU resiues of chain A but not the H atoms
        xyz = db.get('x,y,z',model=0) #chainID='A',resName=['VAL','LEU'],no_name=['H'])

        # put the data back
        db.update('x,y,z',xyz)#,chainID='A',resName=['VAL','LEU'],no_name=['H'])

        # extract atoms
        #atoms = db.get(chainID='A',resName=['VAL','LEU'],no_name=['H'])

        #for at in atoms:
        #   print(at.name,at.x,at.y,at.z)

    def setUp(self):
        self.pdb = 'pdb/test_model.pdb'

if __name__ == '__main__':
    unittest.main()