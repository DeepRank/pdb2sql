import os
import unittest
from pdb2sql.superpose import superpose


class TestSuperpose(unittest.TestCase):
    """Test the superpose functionality"""

    def setUp(self):
        self.pdb1 = 'pdb/1AK4/1AK4_5w.pdb'
        self.pdb2 = 'pdb/1AK4/1AK4_10w.pdb'

    def test_superpose(self):
        """Test superpose()"""
        with self.assertWarns(UserWarning) as ex:
            superpose(self.pdb1, self.pdb2, chainID='A')

    def test_superpose_backbone_error(self):
        """Test superpose() backbone error when specifying `name`"""
        with self.assertWarns(UserWarning) as ex:
            with self.assertRaises(ValueError) as err:
                superpose(self.pdb1, self.pdb2, name='CA')
        err_msg = err.exception.args[0]
        target = "Atom type specified but only_backbone == True"
        self.assertIn(target, err_msg)

if __name__ == '__main__':
    unittest.main()