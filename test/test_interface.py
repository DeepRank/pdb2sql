import unittest
from pathlib import Path
from pdb2sql import interface
from pdb2sql import pdb2sql

from . import pdb_folder


class Test_1_ContactAtoms(unittest.TestCase):
    """Test function get_contact_atoms."""

    def setUp(self):
        self.pdb = Path(pdb_folder, '3CRO.pdb')
        self.db = interface(self.pdb)

    def test_get_contact_atoms_default(self):
        """"verify get_contact_atoms default."""
        contact_atoms = self.db.get_contact_atoms()
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 2)
        self.assertEqual(list(contact_atoms.keys()), ['A', 'B'])
        # in pymol `select natoms, chain A within 8.5 of chain B`
        # to get the number of contact atoms
        self.assertEqual(len(contact_atoms['A']), 341)
        self.assertEqual(len(contact_atoms['B']), 333)

    def test_get_contact_atoms_cutoff(self):
        """"verify get_contact_atoms(cutoff=5.5)"""
        cutoff = 5.5
        contact_atoms = self.db.get_contact_atoms(cutoff=cutoff)
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 2)
        self.assertEqual(list(contact_atoms.keys()), ['A', 'B'])
        self.assertEqual(len(contact_atoms['A']), 185)
        self.assertEqual(len(contact_atoms['B']), 174)

    def test_get_contact_atoms_allchains(self):
        """"verify get_contact_atoms(allchains=True)"""
        contact_atoms = self.db.get_contact_atoms(allchains=True)
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 4)
        self.assertEqual(list(contact_atoms.keys()), ['A', 'B', 'L', 'R'])
        self.assertEqual(len(contact_atoms['A']), 367)
        self.assertEqual(len(contact_atoms['B']), 372)
        self.assertEqual(len(contact_atoms['L']), 314)
        self.assertEqual(len(contact_atoms['R']), 304)

    def test_get_contact_atoms_chain1chain2(self):
        """"verify get_contact_atoms(chain1='L', chain2='R')"""
        contact_atoms = self.db.get_contact_atoms(chain1='L', chain2='R')
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 2)
        self.assertEqual(list(contact_atoms.keys()), ['L', 'R'])
        self.assertEqual(len(contact_atoms['L']), 132)
        self.assertEqual(len(contact_atoms['R']), 132)

    def test_get_contact_atoms_extend2residue(self):
        """"verify get_contact_atoms(extend_to_residue=True)"""
        contact_atoms = self.db.get_contact_atoms(extend_to_residue=True)
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 2)
        self.assertEqual(list(contact_atoms.keys()), ['A', 'B'])
        # in pymol `select natoms, byres(chain A within 8.5 of chain B)`
        # to get the number of contact atoms
        self.assertEqual(len(contact_atoms['A']), 405)
        self.assertEqual(len(contact_atoms['B']), 409)

    def test_get_contact_atoms_onlybackbone_NA(self):
        """"verify get_contact_atoms(extend_to_residue=True) for nuclear
        acids."""
        with self.assertWarns(UserWarning) as ex:
            contact_atoms = self.db.get_contact_atoms(only_backbone_atoms=True)
        self.assertEqual(len(ex.warnings), 1)
        self.assertEqual(ex.warning.args[0],
                         'No contact atoms detected in pdb2sql')
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 2)
        self.assertEqual(list(contact_atoms.keys()), ['A', 'B'])
        self.assertEqual(len(contact_atoms['A']), 0)
        self.assertEqual(len(contact_atoms['B']), 0)

    def test_get_contact_atoms_onlybackbone_protein(self):
        """"verify get_contact_atoms(extend_to_residue=True) for proteins."""
        contact_atoms = self.db.get_contact_atoms(
            only_backbone_atoms=True,
            chain1='L',
            chain2='R'
        )
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 2)
        self.assertEqual(list(contact_atoms.keys()), ['L', 'R'])
        # pymol `select catoms, (chain L and name CA+C+N+O)
        #   within 8.5 of (chain R and name CA+C+N+O)`
        self.assertEqual(len(contact_atoms['L']), 22)
        self.assertEqual(len(contact_atoms['R']), 20)

    def test_get_contact_atoms_exludeH(self):
        """"verify get_contact_atoms(excludeH=True)"""
        pdb = Path(pdb_folder, '3CRO_H.pdb')
        db = interface(pdb)
        contact_atoms = db.get_contact_atoms(excludeH=True)
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 2)
        self.assertEqual(list(contact_atoms.keys()), ['A', 'B'])
        self.assertEqual(len(contact_atoms['A']), 341)
        self.assertEqual(len(contact_atoms['B']), 333)

    def test_get_contact_atoms_contactpairs(self):
        """"verify get_contact_atoms(return_conact_pairs=True)"""
        contact_atoms = self.db.get_contact_atoms(
            return_contact_pairs=True
        )
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 341)

        for i in contact_atoms.keys():
            with self.subTest(i=i):
                self.assertIsInstance(contact_atoms[i], list)
                self.assertNotEqual(len(contact_atoms[i]), 0)
        self.assertEqual(len(contact_atoms[6]), 1)
        self.assertEqual(len(contact_atoms[404]), 19)

    def test_get_contact_atoms_alltrue(self):
        """"verify get_contact_atoms(True)"""
        pdb = Path(pdb_folder, '3CRO_H.pdb')
        db = interface(pdb)
        contact_atoms = db.get_contact_atoms(
            allchains=True,
            extend_to_residue=True,
            only_backbone_atoms=True,
            excludeH=True)
        self.assertIsInstance(contact_atoms, dict)
        self.assertEqual(len(contact_atoms), 4)
        self.assertEqual(list(contact_atoms.keys()), ['A', 'B', 'L', 'R'])
        # pymol `select catoms, name CA+C+N+O and byres((chain L and name CA+C+N+O )
        #  within 8.5 of (chain R and name CA+C+N+O))`
        self.assertEqual(len(contact_atoms['A']), 0)
        self.assertEqual(len(contact_atoms['B']), 0)
        self.assertEqual(len(contact_atoms['L']), 36)
        self.assertEqual(len(contact_atoms['R']), 32)


class Test_2_ContactResidues(unittest.TestCase):
    """test get_contact_residues function."""

    def setUp(self):
        self.pdb = Path(pdb_folder, '3CRO.pdb')
        self.db = interface(self.pdb)

    def test_get_contact_residues_default(self):
        """"verify get_contact_residues default."""
        contact_residues = self.db.get_contact_residues()
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 2)
        self.assertEqual(list(contact_residues.keys()), ['A', 'B'])
        # in pymol:
        #   select natoms, chain A within 8.5 of chain B
        #   stored.nres = set()
        #   iterate (natoms), stored.nres.add((chain, resi, resn))
        #   print(len(stored.nres))
        self.assertEqual(len(contact_residues['A']), 20)
        self.assertEqual(len(contact_residues['B']), 20)

    def test_get_contact_residues_cutoff(self):
        """"verify get_contact_residues(cutoff=5.5)"""
        cutoff = 5.5
        contact_residues = self.db.get_contact_residues(cutoff=cutoff)
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 2)
        self.assertEqual(list(contact_residues.keys()), ['A', 'B'])
        self.assertEqual(len(contact_residues['A']), 20)
        self.assertEqual(len(contact_residues['B']), 20)

    def test_get_contact_residues_allchains(self):
        """"verify get_contact_residues(allchains=True)"""
        contact_residues = self.db.get_contact_residues(allchains=True)
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 4)
        self.assertEqual(list(contact_residues.keys()), ['A', 'B', 'L', 'R'])
        self.assertEqual(len(contact_residues['A']), 20)
        self.assertEqual(len(contact_residues['B']), 20)
        self.assertEqual(len(contact_residues['L']), 47)
        self.assertEqual(len(contact_residues['R']), 48)

    def test_get_contact_residues_chain1chain2(self):
        """"verify get_contact_residues(chain1='L', chain2='R')"""
        contact_residues = self.db.get_contact_residues(chain1='L', chain2='R')
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 2)
        self.assertEqual(list(contact_residues.keys()), ['L', 'R'])
        self.assertEqual(len(contact_residues['L']), 20)
        self.assertEqual(len(contact_residues['R']), 23)

    def test_get_contact_residues_exludeH(self):
        """"verify get_contact_residues(excludeH=True)"""
        pdb = Path(pdb_folder, '3CRO_H.pdb')
        db = interface(pdb)
        contact_residues = db.get_contact_residues(
            allchains=True, excludeH=True)
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 4)
        self.assertEqual(list(contact_residues.keys()), ['A', 'B', 'L', 'R'])
        self.assertEqual(len(contact_residues['A']), 20)
        self.assertEqual(len(contact_residues['B']), 20)
        self.assertEqual(len(contact_residues['L']), 47)
        self.assertEqual(len(contact_residues['R']), 48)

    def test_get_contact_residues_onlybackbone_NA(self):
        """"verify get_contact_residues(only_backbone_atoms=True) for NA."""
        with self.assertWarns(UserWarning) as ex:
            contact_residues = self.db.get_contact_residues(
                only_backbone_atoms=True)
        self.assertEqual(len(ex.warnings), 1)
        self.assertEqual(ex.warning.args[0],
                         'No contact atoms detected in pdb2sql')
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 2)
        self.assertEqual(list(contact_residues.keys()), ['A', 'B'])
        # pymol `select catoms, (chain L and name CA+C+N+O)
        #   within 8.5 of (chain R and name CA+C+N+O)`
        self.assertEqual(len(contact_residues['A']), 0)
        self.assertEqual(len(contact_residues['B']), 0)

    def test_get_contact_residues_onlybackbone_protein(self):
        """"verify get_contact_residues(only_backbone_atoms=True) for
        proteins."""
        contact_residues = self.db.get_contact_residues(
            only_backbone_atoms=True,
            chain1='L',
            chain2='R'
        )
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 2)
        self.assertEqual(list(contact_residues.keys()), ['L', 'R'])
        # pymol `select catoms, (chain L and name CA+C+N+O)
        #   within 8.5 of (chain R and name CA+C+N+O)`
        self.assertEqual(len(contact_residues['L']), 9)
        self.assertEqual(len(contact_residues['R']), 8)

    def test_get_contact_residues_contactpairs(self):
        """"verify get_contact_residues(return_conact_pairs=True)"""
        contact_residues = self.db.get_contact_residues(
            chain1='L', chain2='R', return_contact_pairs=True)
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 20)
        for i in contact_residues.keys():
            with self.subTest(i=i):
                self.assertIsInstance(contact_residues[i], list)
                self.assertNotEqual(len(contact_residues[i]), 0)
        # in pymol:
        # select natoms, (chain R) within 8.5 of (chain L and resi 60)
        self.assertEqual(len(contact_residues[('L', 60, 'GLN')]), 3)

    def test_get_contact_residues_alltrue(self):
        """"verify get_contact_residues(True)"""
        pdb = Path(pdb_folder, '3CRO_H.pdb')
        db = interface(pdb)
        contact_residues = db.get_contact_residues(
            allchains=True, only_backbone_atoms=True, excludeH=True)
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 4)
        self.assertEqual(list(contact_residues.keys()), ['A', 'B', 'L', 'R'])
        self.assertEqual(len(contact_residues['A']), 0)
        self.assertEqual(len(contact_residues['B']), 0)
        self.assertEqual(len(contact_residues['L']), 9)
        self.assertEqual(len(contact_residues['R']), 8)


class Test_3_PDB2SQLInstanceInput(unittest.TestCase):
    """test using pdb2sql instance as input"""

    def setUp(self):
        self.pdb = Path(pdb_folder, '3CRO.pdb')

    def test_get_contact_residues_default(self):
        """"verify get_contact_residues default."""
        pdb_db = pdb2sql(self.pdb)
        self.db = interface(pdb_db)
        contact_residues = self.db.get_contact_residues()
        self.assertIsInstance(contact_residues, dict)
        self.assertEqual(len(contact_residues), 2)
        self.assertEqual(list(contact_residues.keys()), ['A', 'B'])
        self.assertEqual(len(contact_residues['A']), 20)
        self.assertEqual(len(contact_residues['B']), 20)

    def test_database_consistency(self):
        """"verify initilizing interface with updated pdb2sql database"""
        pdb_db = pdb2sql(self.pdb)
        pdb_db.update_column('temp', [99]*10)
        target = pdb_db.get('*')

        self.db = interface(pdb_db)
        result = self.db.get('*')
        self.assertEqual(target, result)

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
