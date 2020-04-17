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

    def test_fetch_nonexist_pdbid(self):
        """Verfify fetch with non-exist PDB ID"""
        pdbid = '1000'
        with self.assertRaisesRegex(ValueError, 'PDB ID not exist'):
            fetch(pdbid)

    def test_fetch_nonexist_pdbfmt(self):
        """Verfify fetch PDB ID that has no pdb format but cif format"""
        pdbid = '6SL9'
        with self.assertRaisesRegex(ValueError,
            'The PDB ID given is only represented in mmCIF format'):
            fetch(pdbid)


if __name__ == "__main__":
    unittest.main()
