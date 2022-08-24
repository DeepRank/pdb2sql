import os
from pathlib import Path
from pdb2sql.StructureSimilarity import StructureSimilarity
import unittest

from . import pdb_folder

class TestSim(unittest.TestCase):
    """Test Similarity calculation."""

    def setUp(self):
        self.decoy = Path(pdb_folder, '1AK4', '1AK4_5w.pdb')
        self.ref = Path(pdb_folder, '1AK4', 'target.pdb')
        self.izone = Path(pdb_folder, '1AK4', 'target.izone')
        self.lzone = Path(pdb_folder, '1AK4', 'target.lzone')
        self.sim = StructureSimilarity(str(self.decoy), str(self.ref))

    #test case for multi-chain pdb files
        self.decoy2 = Path(pdb_folder, 'multi_chain', 'model.pdb')
        self.ref2 = Path(pdb_folder, 'multi_chain', 'ref.pdb')
        self.sim2 = StructureSimilarity(str(self.decoy2), str(self.ref2), enforce_residue_matching = False)
        self.irmsd_multichn = 1.635 #need to double check with ProFit

        # target values are calcualted using scripts from
        # https://github.com/haddocking/BM5-clean
        self.irmsd = 1.135
        self.lrmsd = 6.655
        self.fnat = 0.790698
        self.dockQ = 0.682191
        self.capriClass = 'medium'
        self.nclashes_ref = 4
    ####################################################################
    # test check_residues to see if pdb files match or not
    ####################################################################
    def test_check_residues(self):
        decoy = Path(pdb_folder, '1AK4', '1AK4_5w_nonmatch.pdb')
        with self.assertRaisesRegex(ValueError,
                                'Residue numbering not identical'):
            sim = StructureSimilarity(decoy, self.ref)
            sim.check_residues()

    ####################################################################
    # test i-rmsd
    ####################################################################
    #def test_irmsdfast_default(self):
    #    """verify compute_irmsd_fast()"""
    #    result = self.sim.compute_irmsd_fast()
    #    self.assertEqual(result, self.irmsd)

    def test_irmsdfast_izone(self):
        """verify compute_irmsd_fast(izone='fast.izone)"""
        result = self.sim.compute_irmsd_fast(izone=self.sim.read_zone(self.izone))
        self.assertEqual(result, self.irmsd)

    def test_irmsdfast_method(self):
        """verify compute_irmsd_fast(method='quaternion')"""
        izone = self.sim.compute_izone(chain_pairs=['A:B'])
        result = self.sim.compute_irmsd_fast(izone = izone, method='quaternion')
        self.assertEqual(result, self.irmsd)

    def test_irmsdfast_check(self):
        """verify compute_irmsd_fast(check=False)"""
        result = self.sim.compute_irmsd_fast(izone=self.sim.read_zone(self.izone), check=False)
        self.assertEqual(result, self.irmsd)

    def test_irmsdfast_multichn(self):
        """verify compute_irmsd_fast() on multi-chain pdbs"""
        # MHC: A and B
        # peptide: C
        # TCR: D and E
        izone = self.sim2.compute_izone(chain_pairs=['ABC:DE'])
        result = self.sim2.compute_irmsd_fast(izone=izone , cutoff=10)
        self.assertEqual(result, self.irmsd_multichn)

        izone = self.sim2.compute_izone(chain_pairs=['AB:DE', 'C:DE'])
        result = self.sim2.compute_irmsd_fast(izone=izone , cutoff=10)
        self.assertEqual(result, self.irmsd_multichn)

    def test_irmsdsql_default(self):
        """verify compute_irmsd_pdb2sql()"""
        result = self.sim.compute_irmsd_pdb2sql()
        self.assertEqual(result, self.irmsd)

    def test_irmsdsql_izone(self):
        """verify compute_irmsd_pdb2sql(izoneFL=sql.izone)"""
        result = self.sim.compute_irmsd_pdb2sql(izoneFL=self.izone)
        self.assertEqual(result, self.irmsd)

    def test_irmssql_method(self):
        """verify compute_irmsd_pdb2sql(method='quaternion')"""
        result = self.sim.compute_irmsd_pdb2sql(method='quaternion')
        self.assertEqual(result, self.irmsd)

    def test_irmsdsql_exportpdb(self):
        """verify compute_irmsd_pdb2sql(exportpath='.')"""
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

    def test_lrmsdsql_default(self):
        """verify compute_lrmsd_pdb2sql()"""
        result = self.sim.compute_lrmsd_pdb2sql()
        self.assertEqual(result, self.lrmsd)

    def test_lrmsdsql_method(self):
        """verify compute_lrmsd_pdb2sql(method='quaternion')"""
        result = self.sim.compute_lrmsd_pdb2sql(method='quaternion')
        self.assertEqual(result, self.lrmsd)

    def test_lrmsdsql_exportpdb(self):
        """verify compute_lrmsd_pdb2sql(exportpath='.')"""
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
        result = self.sim.compute_fnat_pdb2sql()
        self.assertEqual(result, self.fnat)

    ####################################################################
    # test dockQ
    ####################################################################
    def test_dockQ_default(self):
        """verify compute_DockQScore()"""
        result = self.sim.compute_DockQScore(self.fnat, self.lrmsd, self.irmsd)
        self.assertEqual(result, self.dockQ)

    ####################################################################
    # test CAPRI
    ####################################################################
    def test_capri_default(self):
        """verify compute_CapriClass()"""
        result = self.sim.compute_CapriClass(self.fnat, self.lrmsd, self.irmsd)
        self.assertEqual(result, self.capriClass)

    def test_capri_dummy(self):
        """verify compute_CapriClass()"""
        fnat = [0.9, 0.8, 0.7, 0.5, 0.3, 0.1]
        lrmsd = [0.8, 2.4, 6.2, 7.5, 12.0, 10.0]
        irmsd = [0.6, 0.8, 1.6, 2.3, 3.1, 3.3]
        targets = ['high', 'high', 'medium',
                   'acceptable', 'acceptable', 'acceptable']
        results = []
        for i, j, k in zip(fnat, lrmsd, irmsd):
            results.append(self.sim.compute_CapriClass(i, j, k))
        self.assertEqual(results, targets)

    ####################################################################
    # test clashes
    ####################################################################
    def test_clashes_default(self):
        """verify compute_clashes()"""
        result = self.sim.compute_clashes(self.ref)
        self.assertEqual(result, self.nclashes_ref)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
