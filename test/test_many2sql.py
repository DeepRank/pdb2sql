import unittest
import os
import numpy as np
from pathlib import Path
from pdb2sql import many2sql
from pdb2sql import pdb2sql
from utils import CaptureOutErr


class TestMany2SQL(unittest.TestCase):

    def setUp(self):
        self.pdb1 = 'pdb/1AK4/1AK4_5w.pdb'
        self.pdb2 = 'pdb/1AK4/1AK4_10w.pdb'

    def test_init_from_files(self):
        """Verify init from path."""
        pdbs = [self.pdb1, self.pdb2]
        many = many2sql(pdbs, tablenames=pdbs)

    def test_init_from_pdb_data(self):
        """Verify init from data."""
        pdbs = [self.pdb1, self.pdb2]
        sqls = [pdb2sql(p) for p in pdbs]
        data = [db.sql2pdb() for db in sqls]
        many = many2sql(data, tablenames=pdbs)

    def test_init_from_sql(self):
        """Verify default sqls."""
        pdbs = [self.pdb1, self.pdb2]
        sqls = [pdb2sql(p) for p in pdbs]
        many = many2sql(sqls, tablenames=pdbs)

    def test_call(self):
        """Test call function."""
        pdbs = [self.pdb1, self.pdb2]
        many = many2sql(pdbs, tablenames=pdbs)
        chainA = many(chainID='A')

    def test_get_all(self):
        """Test get_all function."""
        pdbs = [self.pdb1, self.pdb2]
        many = many2sql(pdbs, tablenames=pdbs)
        data = many.get_all('x,y,z', chainID='A')

    def test_get_intersection(self):
        """Test get_all function."""
        pdbs = [self.pdb1, self.pdb2]
        many = many2sql(pdbs, tablenames=pdbs)
        data = many.get_intersection('x,y,z')

    def test_intersect(self):
        """Test get_all function."""
        pdbs = [self.pdb1, self.pdb2]
        many = many2sql(pdbs, tablenames=pdbs)
        chainA = many.intersect()


if __name__ == '__main__':
    # runner = unittest.TextTestRunner(verbosity=2)
    # unittest.main(testRunner=runner)
    unittest.main()
