import unittest
import os
import numpy as np
from pathlib import Path
from pdb2sql import pdb2sql

from .utils import CaptureOutErr
from . import pdb_folder, test_folder


class TestCreateSQL(unittest.TestCase):

    def setUp(self):
        self.pdbfile = Path(pdb_folder, "3CRO.pdb")
        self.sqlfile = Path(test_folder, "sql_3CRO.db")

        self.pdb_longline = Path(pdb_folder, "dummy_longline.pdb")
        self.pdb_nochainID_segID = Path(pdb_folder, "dummy_blank_chainID_with_segID.pdb")
        self.pdb_nochainID_nosegID = Path(pdb_folder, "dummy_blank_chainID_without_segID.pdb")
        self.pdb_noocc = Path(pdb_folder, "dummy_blank_occupancy.pdb")
        self.pdb_notemp = Path(pdb_folder, "dummy_blank_temperature.pdb")
        self.pdb_noelement = Path(pdb_folder, "dummy_blank_element.pdb")

    def test_init(self):
        """Verify default init."""
        db = pdb2sql(self.pdbfile)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1856)
        last_line = (1859, "O", "", "GLY", "R", 62, "",
                     -32.180, -32.765, 46.907, 1.00, 38.84, "O", 0)
        self.assertEqual(result[-1], last_line)

    def test_pdbfile_Path(self):
        """Verify input as Path instance"""
        p = Path(self.pdbfile)
        db = pdb2sql(p)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1856)
        last_line = (1859, "O", "", "GLY", "R", 62, "",
                     -32.180, -32.765, 46.907, 1.00, 38.84, "O", 0)
        self.assertEqual(result[-1], last_line)

    def test_pdbfile_strings(self):
        """Verify input with a string of pdb content"""
        with open(self.pdbfile) as f:
            pdb = f.read()
        db = pdb2sql(pdb)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1856)
        last_line = (1859, "O", "", "GLY", "R", 62, "",
                     -32.180, -32.765, 46.907, 1.00, 38.84, "O", 0)
        self.assertEqual(result[-1], last_line)

    def test_pdbfile_bytes(self):
        """Verify input with bytes of pdb content"""
        with open(self.pdbfile, 'rb') as f:
            pdb = f.read()
        db = pdb2sql(pdb)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1856)
        last_line = (1859, "O", "", "GLY", "R", 62, "",
                     -32.180, -32.765, 46.907, 1.00, 38.84, "O", 0)
        self.assertEqual(result[-1], last_line)

    def test_pdbfile_list_strings(self):
        """Verify list of strings as input pdb."""
        with open(self.pdbfile) as f:
            pdb = f.readlines()
        db = pdb2sql(pdb)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1856)
        last_line = (1859, "O", "", "GLY", "R", 62, "",
                     -32.180, -32.765, 46.907, 1.00, 38.84, "O", 0)
        self.assertEqual(result[-1], last_line)

    def test_pdbfile_list_bytes(self):
        """Verify list of bytes as input pdb."""

        with open(self.pdbfile) as f:
            pdb = [line.encode() for line in f]
        db = pdb2sql(pdb)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1856)
        last_line = (1859, "O", "", "GLY", "R", 62, "",
                     -32.180, -32.765, 46.907, 1.00, 38.84, "O", 0)
        self.assertEqual(result[-1], last_line)

    def test_pdbfile_array(self):
        """Verify numpy array as input pdb."""

        with open(self.pdbfile) as f:
            pdb = np.array(f.readlines())
        db = pdb2sql(pdb)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1856)
        last_line = (1859, "O", "", "GLY", "R", 62, "",
                     -32.180, -32.765, 46.907, 1.00, 38.84, "O", 0)
        self.assertEqual(result[-1], last_line)

    def test_fix_chainID(self):
        """Verify fix chainID."""
        db = pdb2sql(self.pdbfile, fix_chainID=True)
        db.c.execute('SELECT * FROM ATOM')
        result = db.c.fetchall()
        result_chainIDs = []
        for line in result:
            chainID = line[4]
            if chainID not in result_chainIDs:
                result_chainIDs.append(chainID)
        target = ["A", "B", "C", "D"]
        self.assertEqual(result_chainIDs, target)

    def test_sqlfile_nosave(self):
        """Verify sqlfile without saving."""
        db = pdb2sql(self.pdbfile, sqlfile=self.sqlfile)
        self.assertTrue(os.path.exists(self.sqlfile))
        self.assertTrue(os.path.isfile(self.sqlfile))

    def test_sqlfile_save(self):
        """Verify sqlfile with saving."""
        db = pdb2sql(self.pdbfile, sqlfile=self.sqlfile)
        db._close(rmdb=False)
        self.assertTrue(os.path.exists(self.sqlfile))
        self.assertTrue(os.path.isfile(self.sqlfile))
        os.remove(self.sqlfile)
    
    # The test below is not deemed necessary
    # It was updated since the atom list per residue is now sorted to prevent wrong comparisons
    
    #def test_call(self):
    #    sqldb = pdb2sql(self.pdbfile)
    #    cpy = sqldb(chainID='A')
    #    cpy.c.execute('SELECT * FROM ATOM')
    #    result = cpy.c.fetchall()
    #    last_line = (405, 'C6', '', 'DT', 'A', 20, '', -
    #                 52.817, -4.887, 21.878, 1.0, 2.0, 'C', 0)
    #    self.assertEqual(result[-1], last_line)

    ####################
    # verify pdb format
    ####################
    def test_pdb_longline(self):
        with self.assertRaises(ValueError) as ex:
            db = pdb2sql(self.pdb_longline)
        ex_msg = ex.exception.args[0]
        target = "pdb line is longer than 80:\nATOM      1  O5'  DA A   1     -16.851  -5.543  74.981  1.00 55.62      A    O  1"
        self.assertEqual(ex_msg, target)

    def test_blank_chainID_without_segID(self):
        """Verify blank chainID without segID."""
        with self.assertRaises(ValueError) as ex:
            db = pdb2sql(self.pdb_nochainID_nosegID)
        result = ex.exception.args[0]
        target = 'chainID not found'
        self.assertEqual(result, target)

    def test_blank_chainID_with_segID(self):
        """Verify blank chainID with segID."""
        with self.assertWarns(UserWarning) as ex:
            db = pdb2sql(self.pdb_nochainID_segID)
        self.assertEqual(len(ex.warnings), 24)
        self.assertEqual(ex.warning.args[0],
                         "Missing chainID and set it with segID")
        db.c.execute('SELECT chainID FROM ATOM')
        result = list(set(db.c.fetchall()))
        target = [('L',)]
        self.assertEqual(result, target)

    def test_blank_occ(self):
        """Verify blank occ."""
        db = pdb2sql(self.pdb_noocc)
        db.c.execute('SELECT occ FROM ATOM')
        result = list(set(db.c.fetchall()))
        target = [(1.0,)]
        self.assertEqual(result, target)

    def test_blank_temp(self):
        """Verify blank temperature."""
        db = pdb2sql(self.pdb_notemp)
        db.c.execute('SELECT temp FROM ATOM')
        result = list(set(db.c.fetchall()))
        target = [(10.0,)]
        self.assertEqual(result, target)

    def test_blank_element(self):
        """Verify blank element."""
        db = pdb2sql(self.pdb_noelement)
        db.c.execute('SELECT element FROM ATOM')
        elements = db.c.fetchall()
        result = []
        for i in elements:
            result.append(i[0])
        target = [
            "N",
            "C",
            "C",
            "O",
            "C",
            "C",
            "S",
            "C",
            "N",
            "C",
            "C",
            "O",
            "C",
            "C",
            "C",
            "O",
            "N",
            "N",
            "C",
            "C",
            "O",
            "C",
            "O",
            "C",
            "H",
            "CA"]
        self.assertEqual(result, target)


class TestPrintGetUpdate(unittest.TestCase):
    def setUp(self):
        pdbfile = Path(pdb_folder, "dummy_template.pdb")
        self.db = pdb2sql(pdbfile)

    def test_print(self):
        with CaptureOutErr() as cm:
            self.db.print()
        self.assertEqual(len(cm.stdout), 26)

    def test_get_colnames(self):
        result = self.db.get_colnames()
        self.assertEqual(len(result), 15)
        self.assertEqual(result[0], 'rowID')
        self.assertEqual(result[-1], 'model')

    def test_print_colnames(self):
        with CaptureOutErr() as cm:
            self.db.print_colnames()
        self.assertEqual(len(cm.stdout), 16)
        self.assertEqual(cm.stdout[0], 'Possible column names are:')
        self.assertEqual(cm.stdout[-1], '\tmodel')
        self.assertEqual(cm.stderr, [])

    def test_1_get_columns(self):
        """Verfity get(columns)"""
        result = self.db.get("*")
        self.assertEqual(len(result), 24)
        chainID = set([i[4] for i in result])
        self.assertEqual(chainID, {"L"})

        result = self.db.get(
            "chainID,resSeq,resName,name,x,y,z,model")
        self.assertEqual(len(result), 24)
        chainID = set([i[0] for i in result])
        self.assertEqual(chainID, {"L"})

    def test_1_get_columns_wrongtype(self):
        """Verfity get(columns) wrong column type."""
        with self.assertRaises(TypeError) as ex:
            _ = self.db.get(['x', 'y'])
        ex_msg = ex.exception.args[0]
        target = "argument columns must be str"
        self.assertEqual(ex_msg, target)

    def test_1_get_columns_wrongname(self):
        """Verfity get(columns) wrong column name."""
        with self.assertRaises(ValueError) as ex:
            _ = self.db.get('chain')
        ex_msg = ex.exception.args[0]
        target = "Invalid column name chain"
        self.assertIn(target, ex_msg)

    def test_1_get_kwargs_and(self):
        """Verfity get(**kwargs)"""
        result = self.db.get('*', resName="MET", name=["CA", "CB"])
        self.assertEqual(len(result), 2)
        atomname = set([i[1] for i in result])
        self.assertEqual(atomname, {'CA', 'CB'})
        resname = set([i[3] for i in result])
        self.assertEqual(resname, {'MET'})

    def test_1_get_kwargs_not(self):
        """Verfity get(**kwargs)"""
        result = self.db.get('*', resName="MET", no_name=["CA", "CB"])
        self.assertEqual(len(result), 6)
        atomname = set([i[1] for i in result])
        self.assertNotIn('CA', atomname)
        self.assertNotIn('CB', atomname)
        resname = set([i[3] for i in result])
        self.assertEqual(resname, {'MET'})

    def test_1_get_kwargs_wrongname(self):
        """Verfity get(**kwargs) wrong column name."""
        with self.assertRaises(ValueError) as ex:
            _ = self.db.get('*', no_chain='A')
        ex_msg = ex.exception.args[0]
        target = "Invalid column name chain"
        self.assertIn(target, ex_msg)

    def test_1_get_sql_limit(self):
        """Verfity get() sql varibale limit 999."""
        with CaptureOutErr() as cm:
            with self.assertRaises(ValueError) as ex:
                result = self.db.get(
                    '*', chainID=['L'] * 950, name=["CA"] * 50)
        ex_msg = ex.exception.args[0]
        target = 'Too many SQL variables'
        self.assertEqual(target, ex_msg)
        self.assertIn('SQL Queries can only handle a total of 999 values',
                      cm.stdout[1])

    def test_2_update(self):
        """Verfity update() default."""
        values = np.array([[1., 2., 3.], [4., 5., 6.]])
        self.db.update(
            "x,y,z",
            values=values,
            resName='MET',
            name=[
                'CA',
                'CB'])
        result = self.db.get(
            "resName,name,x,y,z",
            resName='MET',
            name=[
                'CA',
                'CB'])
        target = [['MET', 'CA', 1.0, 2.0, 3.0],
                  ['MET', 'CB', 4.0, 5.0, 6.0]]
        self.assertEqual(result, target)

    def test_2_update_column_wrongtype(self):
        """Verfity update() wrong column type."""
        values = np.array([[1., 2., 3.], [4., 5., 6.]])
        with self.assertRaises(TypeError) as ex:
            self.db.update(['x', 'y', 'z'], values=values,
                           resName='MET', name=['CA', 'CB'])
        ex_msg = ex.exception.args[0]
        target = "argument columns must be str"
        self.assertEqual(ex_msg, target)

    def test_2_update_column_wrongname(self):
        """Verfity update() wrong column name."""
        values = np.array([[1., 2., 3.], [4., 5., 6.]])
        with self.assertRaises(ValueError) as ex:
            self.db.update('xyz', values=values,
                           resName='MET', name=['CA', 'CB'])
        ex_msg = ex.exception.args[0]
        target = "Invalid column name xyz"
        self.assertIn(target, ex_msg)

    def test_2_update_columns_values_mismatch(self):
        """Verfity update() columns and values mismatch."""
        values = np.array([[1., 2., 3.], [4., 5., 6.]])
        with self.assertRaises(ValueError) as ex:
            self.db.update('x,y', values=values,
                           resName='MET', name=['CA', 'CB'])
        ex_msg = ex.exception.args[0]
        target = 'Number of cloumns does not match between argument columns and values'
        self.assertIn(target, ex_msg)

    def test_2_update_rows_values_mismatch(self):
        """Verfity update() rows and values mismatch."""
        values = np.array([[1., 2., 3.], [4., 5., 6.]])
        with self.assertRaises(ValueError) as ex:
            self.db.update('x,y,z', values=values,
                           resName='MET', name=['CA'])
        ex_msg = ex.exception.args[0]
        target = 'Number of data values incompatible with the given conditions'
        self.assertIn(target, ex_msg)

    def test_3_update_cloumn(self):
        """Verfity update_column() default."""
        values = list(range(24))
        self.db.update_column("x", values=values)
        result = self.db.get("x")
        target = values
        self.assertEqual(result, target)

    def test_3_update_cloumn_index(self):
        """Verfity update_column() default."""
        values = [3., 2., 1.]
        # Values sorted by: chainID, resSeq, name
        values_sorted = [1., 2., 3.]
        self.db.update_column(
            "x", values=values, index=list(range(3)))
        result = self.db.get("x", resName='MET',
                             name=['N', 'CA', 'C'])
        # The returned values will be for ['C', 'CA', 'N']
        target = values_sorted
        self.assertEqual(result, target)

    def test_4_add_cloumn(self):
        """Verfity add_column() default."""
        self.db.add_column("id")
        result = self.db.get("id")
        self.assertEqual(len(result), 24)
        self.assertEqual(set(result), {0.})

    def test_4_add_cloumn_str(self):
        """Verfity add_column() add string."""
        values = 'C'
        self.db.add_column("id", coltype='str', value=values)
        result = self.db.get("id")
        self.assertEqual(len(result), 24)
        self.assertEqual(set(result), {'C'})


if __name__ == '__main__':
    # runner = unittest.TextTestRunner(verbosity=2)
    # unittest.main(testRunner=runner)
    unittest.main()
