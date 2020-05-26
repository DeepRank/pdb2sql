import unittest
import numpy as np
from pathlib import Path
from pdb2sql import many2sql
from pdb2sql import pdb2sql

from .utils import CaptureOutErr
from . import pdb_folder


class TestMany2SQL(unittest.TestCase):

    def setUp(self):
        pdb1 = Path(pdb_folder, '1AK4', '1AK4_5w.pdb')
        pdb2 = Path(pdb_folder, '1AK4', '1AK4_10w.pdb')
        self.pdbs = [pdb1, pdb2]
        self.tablenames = [str(pdb1), str(pdb2)]

    def test_init_from_files(self):
        """Verify init from path."""
        many = many2sql(self.pdbs, tablenames=self.tablenames)

    def test_init_from_pdb_data(self):
        """Verify init from data."""
        sqls = [pdb2sql(p) for p in self.pdbs]
        data = [db.sql2pdb() for db in sqls]
        many = many2sql(data, tablenames=self.tablenames)

    def test_init_from_sql(self):
        """Verify default sqls."""
        sqls = [pdb2sql(p) for p in self.pdbs]
        many = many2sql(sqls, tablenames=self.tablenames)

    def test_call(self):
        """Test call function."""
        many = many2sql(self.pdbs, tablenames=self.tablenames)
        chainA = many(chainID='A')

    def test_get_all(self):
        """Test get_all function."""
        many = many2sql(self.pdbs, tablenames=self.tablenames)
        data = many.get_all('x,y,z', chainID='A')

    def test_get_intersection(self):
        """Test get_all function."""
        many = many2sql(self.pdbs, tablenames=self.tablenames)
        data = many.get_intersection('x,y,z')

    def test_intersect(self):
        """Test get_all function."""
        many = many2sql(self.pdbs, tablenames=self.tablenames)
        chainA = many.intersect()


if __name__ == '__main__':
    unittest.main()