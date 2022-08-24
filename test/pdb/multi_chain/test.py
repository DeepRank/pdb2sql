import pdb2sql

interface = pdb2sql.interface('model.pdb')
xyz = interface.get_contact_atoms(chain1='AB',chain2='D')
print(xyz)