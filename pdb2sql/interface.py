
import numpy as np
import itertools
from .pdb2sqlcore import pdb2sql

#from pdb2sqlAlchemy import pdb2sql_alchemy as pdb2sql

'''
Class that allows to analyze itnerface between two chains
of a given complex

Works with pdb2sql AND pdb2sql_alchemy
Some methods an be made simpler using alcemy specific queries

We could also returns a series of OBJECTS instead of the raw data

'''


class interface(pdb2sql):

    def __init__(self,pdb):

        super().__init__(pdb,no_extra=True)
        super()._create_sql()
        self.backbone_type  = ['CA','C','N','O']

    ############################################################################
    #
    # get the contact atoms
    #
    #############################################################################

    def get_contact_atoms(self,cutoff=8.5,allchains=False,chain1='A',chain2='B',
                          extend_to_residue=False,only_backbone_atoms=False,
                          excludeH=False,return_only_backbone_atoms=False,return_contact_pairs=False):

        if allchains:
            chainIDs = super().get_chains()
        else:
            chainIDs = [chain1,chain2]
        nchains = len(chainIDs)

        xyz = dict()
        index = dict()
        resName = dict()
        atName = dict()

        for chain in chainIDs:

            print(chain)
            data = np.array(super().get('x,y,z,rowID,resName,name',chainID=chain))
            print(data)
            xyz[chain] = data[:,:3].astype(float)
            index[chain] = data[:,3].astype(int)
            resName[chain] = data[:,-2]
            atName[chain] = data[:,-1]

        # # xyz of the chains
        # xyz1 = np.array(super().get('x,y,z',chainID=chain1))
        # xyz2 = np.array(super().get('x,y,z',chainID=chain2))

        # # index of b
        # index1 = super().get('rowID',chainID=chain1)
        # index2 = super().get('rowID',chainID=chain2)

        # # resName of the chains
        # resName1 = np.array(super().get('resName',chainID=chain1))
        # resName2 = np.array(super().get('resName',chainID=chain2))

        # # atomnames of the chains
        # atName1 = np.array(super().get('name',chainID=chain1))
        # atName2 = np.array(super().get('name',chainID=chain2))


        # loop through the first chain
        # TO DO : loop through the smallest chain instead ...
        index_contact_1,index_contact_2 = [],[]
        index_contact_pairs = {}

        index_contact = dict()
        index_contact_pairs = {}

        for chain1,chain2 in itertools.combinations(chainIDs,2):

            xyz1 = xyz[chain1]
            xyz2 = xyz[chain2]

            atName1 = atName[chain1]
            atName2 = atName[chain2]


            if chain1 not in index_contact:
                index_contact[chain1] = []

            if chain2 not in index_contact:
                index_contact[chain2] = []

            for i,x0 in enumerate(xyz1):

                # compute the contact atoms
                contacts = np.where(np.sqrt(np.sum((xyz2-x0)**2,1)) <= cutoff )[0]

                # exclude the H if required
                if excludeH and atName1[i][0] == 'H':
                    continue

                if len(contacts)>0 and any([not only_backbone_atoms, atName1[i] in self.backbone_type]):

                    # the contact atoms
                    index_contact[chain1] += [index[chain1][i]]
                    index_contact[chain2] += [index[chain2][k] for k in contacts if ( any( [atName2[k] in self.backbone_type,  not only_backbone_atoms]) and not (excludeH and atName2[k][0]=='H') ) ]

                    # the pairs
                    pairs = [index[chain2][k] for k in contacts if any( [atName2[k] in self.backbone_type,  not only_backbone_atoms] ) and not (excludeH and atName2[k][0]=='H') ]
                    if len(pairs) > 0:
                        index_contact_pairs[index[chain1][i]] = pairs

        # get uniques
        for chain in chainIDs:
            index_contact[chain] = sorted(set(index_contact[chain]))


        # if no atoms were found
        if len(index_contact_pairs)==0:
            print('Warning : No contact atoms detected in pdb2sql')

        # extend the list to entire residue
        if extend_to_residue:
            for chain in chainIDs:
                index_contact[chain] = self._extend_contact_to_residue(index_contact_1,only_backbone_atoms)
            #index_contact_1,index_contact_2 = self._extend_contact_to_residue(index_contact_1,index_contact_2,only_backbone_atoms)


        # filter only the backbone atoms
        if return_only_backbone_atoms and not only_backbone_atoms:

            # get all the names
            # there are better ways to do that !
            atNames = np.array(super().get('name'))

            # change the index_contacts
            for chain in chainIDs:
                index_contact[chain] = [  ind for ind in index_contact[chain] if atNames[ind] in self.backbone_type ]


            # change the contact pairs
            tmp_dict = {}
            for ind1,ind2_list in index_contact_pairs.items():

                if atNames[ind1] in self.backbone_type:
                    tmp_dict[ind1] = [ind2 for ind2 in ind2_list if atNames[ind2] in self.backbone_type]

            index_contact_pairs = tmp_dict

        # not sure that's the best way of dealing with that
        if return_contact_pairs:
            return index_contact_pairs
        else:
            return index_contact

    # extend the contact atoms to the residue
    def _extend_contact_to_residue(self,index1,only_backbone_atoms):

        # extract the data
        dataA = super().get('chainID,resName,resSeq',rowID=index1)
        #dataB = super().get('chainID,resName,resSeq',rowID=index2)

        # create tuple cause we want to hash through it
        dataA = list(map(lambda x: tuple(x),dataA))
        #dataB = list(map(lambda x: tuple(x),dataB))

        # extract uniques
        resA = list(set(dataA))
        #resB = list(set(dataB))

        # init the list
        index_contact_A = []

        # contact of chain A
        for resdata in resA:
            chainID,resName,resSeq = resdata

            if only_backbone_atoms:
                index = super().get('rowID',chainID=chainID,resName=resName,resSeq=resSeq)
                name = super().get('name',chainID=chainID,resName=resName,resSeq=resSeq)
                index_contact_A += [ ind for ind,n in zip(index,name) if n in self.backbone_type ]
            else:
                index_contact_A += super().get('rowID',chainID=chainID,resName=resName,resSeq=resSeq)

        # # contact of chain B
        # for resdata in resB:
        #     chainID,resName,resSeq = resdata
        #     if only_backbone_atoms:
        #         index = self.get('rowID',chainID=chainID,resName=resName,resSeq=resSeq)
        #         name = self.get('name',chainID=chainID,resName=resName,resSeq=resSeq)
        #         index_contact_B += [ ind for ind,n in zip(index,name) if n in self.backbone_type ]
        #     else:
        #         index_contact_B += super().get('rowID',chainID=chainID,resName=resName,resSeq=resSeq)

        # make sure that we don't have double (maybe optional)
        index_contact_A = sorted(set(index_contact_A))
        #index_contact_B = sorted(set(index_contact_B))

        return index_contact_A #,index_contact_B


    # get the contact residue
    def get_contact_residues(self,cutoff=8.5,allchains=False,chain1='A',chain2='B',excludeH=False,
                            only_backbone_atoms=False,return_contact_pairs=False):

        # get the contact atoms
        if return_contact_pairs:

            # declare the dict
            residue_contact_pairs = {}

            # get the contact atom pairs
            atom_pairs = self.get_contact_atoms(cutoff=cutoff,allchains=allchains,chain1=chain1,chain2=chain2,
                                                only_backbone_atoms=only_backbone_atoms,
                                                excludeH=excludeH,
                                                return_contact_pairs=True)

            # loop over the atom pair dict
            for iat1,atoms2 in atom_pairs.items():

                # get the res info of the current atom
                data1 = tuple(super().get('chainID,resSeq,resName',rowID=[iat1])[0])

                # create a new entry in the dict if necessary
                if data1 not in residue_contact_pairs:
                    residue_contact_pairs[data1] = set()

                # get the res info of the atom in the other chain
                data2 = super().get('chainID,resSeq,resName',rowID=atoms2)

                # store that in the dict without double
                for resData in data2:
                    residue_contact_pairs[data1].add(tuple(resData))

            for resData in residue_contact_pairs.keys():
                residue_contact_pairs[resData] = sorted(residue_contact_pairs[resData])

            return residue_contact_pairs

        else:

            # get the contact atoms
            contact_atoms = self.get_contact_atoms(cutoff=cutoff,allchains=allchains,chain1=chain1,chain2=chain2,
                                                   return_contact_pairs=False)

            # get the residue info
            data = dict()
            residue_contact = dict()
            for chain in contact_atoms.keys():
                print(chain)
                print(type(contact_atoms[chain]))
                print(contact_atoms[chain])
                data[chain] = super().get('chainID,resSeq,resName',rowID=list(contact_atoms[chain]))
                print(data[chain])
                residue_contact[chain] = sorted(set([tuple(resData) for resData in data[chain]]))


            return residue_contact


