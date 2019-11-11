import os
from pdb2sql.StructureSimilarity import StructureSimilarity
import unittest


class TestSim(unittest.TestCase):
    """Test Similarity calculation."""

    def setUp(self):
        self.decoy = 'pdb/1AK4/1AK4_5w.pdb'
        self.ref = 'pdb/1AK4/target.pdb'
        self.izone = 'pdb/1AK4/target.izone'
        self.lzone = 'pdb/1AK4/target.lzone'
        self.sim = StructureSimilarity(self.decoy, self.ref)
        # target values are calcualted using scripts from
        # https://github.com/haddocking/BM5-clean
        self.irmsd = 1.135
        self.lrmsd = 6.655
        self.fnat = 0.790698
        self.dockQ = 0.682191

    ####################################################################
    # test i-rmsd
    ####################################################################
    def test_irmsdfast_default(self):
        """verify compute_irmsd_fast()"""
        result = self.sim.compute_irmsd_fast()
        self.assertEqual(result, self.irmsd)

    def test_irmsdfast_izone(self):
        """verify compute_irmsd_fast(izone='fast.izone)"""
        result = self.sim.compute_irmsd_fast(izone=self.izone)
        self.assertEqual(result, self.irmsd)

    def test_irmsdfast_method(self):
        """verify compute_irmsd_fast(method='quaternion')"""
        result = self.sim.compute_irmsd_fast(method='quaternion')
        self.assertEqual(result, self.irmsd)

    def test_irmsdfast_check(self):
        """verify compute_irmsd_fast(check=False)"""
        result = self.sim.compute_irmsd_fast(check=False)
        self.assertEqual(result, self.irmsd)

    def test_irmsdsql_default(self):
        """verify compute_irmsd_pdb2sql()"""
        with self.assertWarns(UserWarning) as ex:
            result = self.sim.compute_irmsd_pdb2sql()
        self.assertEqual(result, self.irmsd)

    def test_irmsdsql_izone(self):
        """verify compute_irmsd_pdb2sql(izone='sql.izone)"""
        with self.assertWarns(UserWarning) as ex:
            result = self.sim.compute_irmsd_pdb2sql(izone=self.izone)
        self.assertEqual(result, self.irmsd)

    def test_irmssql_method(self):
        """verify compute_irmsd_pdb2sql(method='quaternion')"""
        with self.assertWarnsRegex(UserWarning, 'Missing element'):
            result = self.sim.compute_irmsd_pdb2sql(method='quaternion')
        self.assertEqual(result, self.irmsd)

    def test_irmsdsql_exportpdb(self):
        """verify compute_irmsd_pdb2sql(exportpath='.')"""
        with self.assertWarns(UserWarning) as ex:
            result = self.sim.compute_irmsd_pdb2sql(exportpath='.')
        self.assertEqual(result, self.irmsd)
        self.assertTrue(os.path.isfile('./irmsd_ref.pdb'))
        self.assertTrue(os.path.isfile('./irmsd_decoy.pdb'))
        self.assertTrue(os.path.getsize('./irmsd_ref.pdb') > 0)
        self.assertTrue(os.path.getsize('./irmsd_decoy.pdb') > 0)
        os.remove('./irmsd_ref.pdb')
        os.remove('./irmsd_decoy.pdb')

    ####################################################################
    # test l-rmsd
    ####################################################################
    def test_lrmsdfast_default(self):
        """verify compute_lrmsd_fast()"""
        result = self.sim.compute_lrmsd_fast()
        self.assertEqual(result, self.lrmsd)

    def test_lrmsdfast_lzone(self):
        """verify compute_lrmsd_fast(lzone='fast.lzone)"""
        result = self.sim.compute_lrmsd_fast(lzone=self.lzone)
        self.assertEqual(result, self.lrmsd)

    def test_lrmsdfast_method(self):
        """verify compute_lrmsd_fast(method='quaternion')"""
        result = self.sim.compute_lrmsd_fast(method='quaternion')
        self.assertEqual(result, self.lrmsd)

    def test_lrmsdfast_check(self):
        """verify compute_lrmsd_fast(check=False)"""
        with self.assertRaisesRegex(ValueError,
            'operands could not be broadcast') as ex:
            result = self.sim.compute_lrmsd_fast(check=False)

    def test_lrmsdsql_default(self):
        """verify compute_lrmsd_pdb2sql()"""
        with self.assertWarns(UserWarning) as ex:
            result = self.sim.compute_lrmsd_pdb2sql()
        self.assertEqual(result, self.lrmsd)

    def test_lrmsdsql_method(self):
        """verify compute_lrmsd_pdb2sql(method='quaternion')"""
        with self.assertWarns(UserWarning) as ex:
            result = self.sim.compute_lrmsd_pdb2sql(method='quaternion')
        self.assertEqual(result, self.lrmsd)

    def test_lrmsdsql_exportpdb(self):
        """verify compute_lrmsd_pdb2sql(exportpath='.')"""
        with self.assertWarnsRegex(UserWarning, 'Missing element'):
            result = self.sim.compute_lrmsd_pdb2sql(exportpath='.')
        self.assertEqual(result, self.lrmsd)
        self.assertTrue(os.path.isfile('./lrmsd_ref.pdb'))
        self.assertTrue(os.path.isfile('./lrmsd_decoy.pdb'))
        self.assertTrue(os.path.getsize('./lrmsd_ref.pdb') > 0)
        self.assertTrue(os.path.getsize('./lrmsd_decoy.pdb') > 0)
        os.remove('./lrmsd_ref.pdb')
        os.remove('./lrmsd_decoy.pdb')

    ####################################################################
    # test FNAT
    ####################################################################
    def test_fnatfast_default(self):
        """verify compute_fnat_fast()"""
        result = self.sim.compute_fnat_fast()
        self.assertEqual(result, self.fnat)

    def test_fnatsql_default(self):
        """verify compute_fnat_pdb2sql()"""
        with self.assertWarnsRegex(UserWarning, 'Missing element'):
            result = self.sim.compute_fnat_pdb2sql()
        self.assertEqual(result, self.fnat)

    ####################################################################
    # test dockQ
    ####################################################################
    def test_dockQ_default(self):
        """verify compute_DockQScore()"""
        result = self.sim.compute_DockQScore(self.fnat, self.lrmsd, self.irmsd)
        self.assertEqual(result, self.dockQ)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
