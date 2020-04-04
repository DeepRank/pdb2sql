import sys
from pdb2sql import align
from pdb2sql import align_interface

pdb = '../test/1AK4/decoys/1AK4_cm-it0_745.pdb'

#- example 1
align(pdb) # align PC1 to axis x

#- example 2
align(pdb, axis = 'z', export = True)

#- example 3
selection = {'chainID':['A'], 'resSeq':['30', '144'], 'name' : ['CA']}
align(pdb,  export = True, **selection)
sys.exit()

# other examples for selection
selection = {'no_chainID':['A'], 'no_name' : ['CA','C', 'O', 'N'], 'no_resName' : ['ALA', 'TRP']}

#- example 4
pdb1 = '../test/1AK4/decoys/1AK4_cm-it0_745.pdb'
pdb2 = '../test/1AK4/decoys/1AK4_cm-itw_238w.pdb'
align_interface(pdb1,  export = True)
align_interface(pdb2,  export = True)

#- example 4


