from pdb2sql.interface import interface
import unittest

class TestSim(unittest.TestCase):
	"""Test Interface calculation."""

	def test_interface(self):


		db = interface(self.pdb)

		#db.get_contact_atoms()
		db.get_contact_residues(chain1='A',chain2='C',allchains=True)

	def setUp(self):
		self.pdb = 'pdb/target136-scoring_0506_conv.pdb'

if __name__ == '__main__':
    unittest.main()