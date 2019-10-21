import time
from pdb2sql.StructureSimilarity import StructureSimilarity
import unittest


class TestSim(unittest.TestCase):
    """Test Similarity calculation."""

    def test_sim(self):

        sim = StructureSimilarity(self.decoy, self.ref)

        # ----------------------------------------------------------------------

        t0 = time.time()
        irmsd_fast = sim.compute_irmsd_fast(method='svd', izone='1AK4.izone')
        t1 = time.time() - t0
        print('\nIRMSD TIME FAST %f in %f sec' % (irmsd_fast, t1))

        t0 = time.time()
        irmsd = sim.compute_irmsd_pdb2sql(method='svd', izone='1AK4.izone')
        t1 = time.time() - t0
        print('IRMSD TIME SQL %f in %f sec' % (irmsd, t1))

        # ----------------------------------------------------------------------

        t0 = time.time()
        lrmsd_fast = sim.compute_lrmsd_fast(
            method='svd', lzone='1AK4.lzone', check=True)
        t1 = time.time() - t0
        print('\nLRMSD TIME FAST %f in %f sec' % (lrmsd_fast, t1))

        t0 = time.time()
        lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None, method='svd')
        t1 = time.time() - t0
        print('LRMSD TIME SQL %f in %f sec' % (lrmsd, t1))

        # ----------------------------------------------------------------------

        t0 = time.time()
        Fnat = sim.compute_fnat_pdb2sql()
        t1 = time.time() - t0
        print('\nFNAT TIME SQL %f in %f sec' % (Fnat, t1))

        t0 = time.time()
        Fnat_fast = sim.compute_fnat_fast(ref_pairs='1AK4.ref_pairs')
        t1 = time.time() - t0
        print('LRMSD TIME FAST %f in %f sec' % (Fnat_fast, t1))

        # ----------------------------------------------------------------------

        dockQ = sim.compute_DockQScore(Fnat_fast, lrmsd_fast, irmsd_fast)
        print('\nDockQ  %f' % dockQ)

    def setUp(self):
        self.decoy = 'pdb/1AK4_324w.pdb'
        self.ref = 'pdb/1AK4.pdb'


if __name__ == '__main__':
    unittest.main()
