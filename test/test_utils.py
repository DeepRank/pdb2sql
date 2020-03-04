import os
import unittest
from pdb2sql import fetch

class TestTools(unittest.TestCase):

    def test_fetch(self):
        """Verfify fetch with valid pdb"""
        pdbid = '1cbh'
        fetch(pdbid)
        self.assertTrue(os.path.isfile('./1cbh.pdb'))
        os.remove('./1cbh.pdb')

    def test_fetch_invalid_pdb(self):
        """Verfify fetch with invalid pdb"""
        pdbids = ['a', 'ab', 'abc', '1cbha', '1cb*', '-1cb']
        for pdbid in pdbids:
            with self.subTest(pdbid=pdbid):
                with self.assertRaisesRegex(ValueError, 'Invalid PDB ID'):
                    fetch(pdbid)

    def test_fetch_nonexist_pdb(self):
        """Verfify fetch with non-exist pdb"""
        pdbid = '1000'
        with self.assertRaisesRegex(ValueError, 'PDB not exist'):
            fetch(pdbid)


if __name__ == "__main__":
    unittest.main()
