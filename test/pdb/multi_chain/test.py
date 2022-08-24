import pdb2sql
from pdb2sql.StructureSimilarity import StructureSimilarity

interface = pdb2sql.interface('model.pdb')
xyz = interface.get_contact_atoms(chain1='AB',chain2='D')
print(xyz)

ss = StructureSimilarity('model.pdb','ref.pdb',  enforce_residue_matching = False)
ss.compute_irmsd_fast(chain_pairs='AB:D')