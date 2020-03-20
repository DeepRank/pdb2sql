import os 
import unittest
from pdb2sql.superpose import superpose

class TestSuperpose(unittest.TestCase):
    """Test the superpose functionality"""

    def setUp(self):

        self.pdb1 = 'pdb/1AK4/1AK4_5w.pdb'
        self.pdb2 = 'pdb/1AK4/1AK4_10w.pdb'

    def test_superpose(self):
        """test the routine."""
        superpose(self.pdb1, self.pdb2, chainID='A')

    @unittest.expectedFailure    
    def test_fails(self):
        superpose(self.pdb1, self.pdb2, name='CA')