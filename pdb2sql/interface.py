import numpy as np
import itertools
import warnings
from .pdb2sqlcore import pdb2sql


class interface(pdb2sql):

    def __init__(self, pdb, **kwargs):
        """Create an independent SQL database for interface object.

        Args:
            pdb(str, list, ndarray, pdb2sql): pdb file or data, or pdb2sql object.
                If pdb2sql object is used, all changes in the database of pdb2sql
                object before initializing the interface instance will be used in the
                new sql database of the interface instance; afterwards, two databses
                will be independent from each other.

        Examples:
            >>> from pdb2sql import pdb2sql
            >>> from pdb2sql import interface
            >>> # use pdb2sql object as input
            >>> pdb_db = pdb2sql('3CRO.pdb')
            >>> interface_db1 = interface(pdb_db)
            >>> # use pdb file as input
            >>> interface_db2 = interface('3CRO.pdb')
        """
        if isinstance(pdb, pdb2sql):
            pdb._commit()
            pdb = pdb.sql2pdb()
        super().__init__(pdb, **kwargs)

    def __repr__(self):
        return f'{self.__module__}.{self.__class__.__name__} object'
    ##########################################################################
    #
    # get the contact atoms
    #
    ##########################################################################

    def get_contact_atoms(
            self,
            cutoff=8.5,
            allchains=False,
            chain1='A',
            chain2='B',
            extend_to_residue=False,
            only_backbone_atoms=False,
            excludeH=False,
            return_contact_pairs=False):
        """Get rowIDs of contact atoms.

        Args:
            cutoff (float): distance cutoff for calculating contact.
                Defaults to 8.5.
            allchains (bool): calculate contacts for all chains or not.
                 Defaults to False.
            chain1 (str): first chain ID. Defaults to 'A'.
                Used when 'allchains' is False.
            chain2 (str): second chain ID. Defaults to 'B'.
                Used when 'allchains' is False.
            extend_to_residue (bool): get all atoms of the residues containing
                at least one contact atom. Defaults to False.
            only_backbone_atoms (bool): only use backbone atoms to
                calculate contact or not. Defaults to False.
            excludeH (bool): Exculde hydrogen atoms for contact
                calculation or not. Defaults to False.
            return_contact_pairs (bool): if return atomic contact pairs
                or not. Defaults to False.

        Returns:
            dict: rowID of contact atoms or rowID of contact atom pairs
        """
        if allchains:
            chainIDs = self.get_chains()
        else:
            chainIDs = [chain1, chain2]

        chains = self.get_chains()
        for c in list(chain1)+list(chain2):
            if c not in chains:
                raise ValueError(
                    'chain %s not found in the structure' % c)

        xyz = dict()
        index = dict()
        resName = dict()
        atName = dict()

        for chain in chainIDs:

            data = np.array(
                self.get('x,y,z,rowID,resName,name', chainID=list(chain)))
            xyz[chain] = data[:, :3].astype(float)
            index[chain] = data[:, 3].astype(int)
            resName[chain] = data[:, -2]
            atName[chain] = data[:, -1]

        # loop through the first chain
        # TODO : loop through the smallest chain instead ...
        #index_contact_1,index_contact_2 = [],[]
        #index_contact_pairs = {}

        index_contact = dict()
        index_contact_pairs = dict()

        for chain1, chain2 in itertools.combinations(chainIDs, 2):

            xyz1 = xyz[chain1]
            xyz2 = xyz[chain2]

            atName1 = atName[chain1]
            atName2 = atName[chain2]

            if chain1 not in index_contact:
                index_contact[chain1] = []

            if chain2 not in index_contact:
                index_contact[chain2] = []

            for i, x0 in enumerate(xyz1):

                # compute the contact atoms
                contacts = np.where(
                    np.sqrt(np.sum((xyz2 - x0)**2, 1)) <= cutoff)[0]

                # exclude the H if required
                if excludeH and atName1[i][0] == 'H':
                    continue

                if len(contacts) > 0 and any(
                        [not only_backbone_atoms, atName1[i] in self.backbone_atoms]):

                    pairs = [
                        index[chain2][k] for k in contacts if any(
                            [
                                atName2[k] in self.backbone_atoms,
                                not only_backbone_atoms]) and not (
                            excludeH and atName2[k][0] == 'H')]
                    if len(pairs) > 0:
                        index_contact_pairs[index[chain1][i]] = pairs
                        index_contact[chain1] += [index[chain1][i]]
                        index_contact[chain2] += pairs

        # if no atoms were found
        if len(index_contact_pairs) == 0:
            warnings.warn('No contact atoms detected in pdb2sql')

        # get uniques
        for chain in chainIDs:
            index_contact[chain] = sorted(set(index_contact[chain]))

        # extend the list to entire residue
        if extend_to_residue:
            for chain in chainIDs:
                index_contact[chain] = self._extend_contact_to_residue(
                    index_contact[chain], only_backbone_atoms)

        # not sure that's the best way of dealing with that
        # TODO split to two functions get_contact_atoms and
        # get_contact_atom_pairs
        if return_contact_pairs:
            return index_contact_pairs
        else:
            return index_contact

    # extend the contact atoms to the residue
    def _extend_contact_to_residue(self, index1, only_backbone_atoms):

        # extract the data
        dataA = self.get('chainID,resName,resSeq', rowID=index1)
        #dataB = self.get('chainID,resName,resSeq',rowID=index2)

        # create tuple cause we want to hash through it
        dataA = list(map(lambda x: tuple(x), dataA))
        #dataB = list(map(lambda x: tuple(x),dataB))

        # extract uniques
        resA = list(set(dataA))
        #resB = list(set(dataB))

        # init the list
        index_contact_A = []

        # contact of chain A
        for resdata in resA:
            chainID, resName, resSeq = resdata

            if only_backbone_atoms:
                index = self.get(
                    'rowID',
                    chainID=chainID,
                    resName=resName,
                    resSeq=resSeq)
                name = self.get(
                    'name',
                    chainID=chainID,
                    resName=resName,
                    resSeq=resSeq)
                index_contact_A += [ind for ind,
                                    n in zip(index,
                                             name) if n in self.backbone_atoms]
            else:
                index_contact_A += self.get('rowID',
                                            chainID=chainID,
                                            resName=resName,
                                            resSeq=resSeq)

        # make sure that we don't have double (maybe optional)
        index_contact_A = sorted(set(index_contact_A))

        return index_contact_A

    # get the contact residue
    def get_contact_residues(
            self,
            cutoff=8.5,
            allchains=False,
            chain1='A',
            chain2='B',
            excludeH=False,
            only_backbone_atoms=False,
            return_contact_pairs=False):
        """Get contact residues represented with (chain,resSeq, resname).

        Args:
            cutoff (float): distance cutoff for contact calculation
                Defaults to 8.5.
            allchains (bool): calculate contacts for all chains or not.
                 Defaults to False.
            chain1 (str): first chain ID. Defaults to 'A'.
            chain2 (str): second chain ID. Defaults to 'B'.
            excludeH (bool): Exculde hydrogen atoms for contact
                calculation or not. Defaults to False.
            only_backbone_atoms (bool): only use backbone atoms to
                calculate contact or not. Defaults to False.
            return_contact_pairs (bool): if return residue contact pairs
                or not. Defaults to False.

        Returns:
            dict: (chain,resSeq,resName) of contact residues or
                contact residue pairs.
        """
        # TODO split this func to two functions
        # TODO get_contact_residues and get_contact_residue_pairs

        # get the contact atoms
        if return_contact_pairs:

            # declare the dict
            residue_contact_pairs = {}

            # get the contact atom pairs
            atom_pairs = self.get_contact_atoms(
                cutoff=cutoff,
                allchains=allchains,
                chain1=chain1,
                chain2=chain2,
                only_backbone_atoms=only_backbone_atoms,
                excludeH=excludeH,
                return_contact_pairs=True)

            # loop over the atom pair dict
            for iat1, atoms2 in atom_pairs.items():

                # get the res info of the current atom
                data1 = tuple(
                    self.get(
                        'chainID,resSeq,resName',
                        rowID=[iat1])[0])

                # create a new entry in the dict if necessary
                if data1 not in residue_contact_pairs:
                    residue_contact_pairs[data1] = set()

                # get the res info of the atom in the other chain
                data2 = self.get(
                    'chainID,resSeq,resName', rowID=atoms2)

                # store that in the dict without double
                for resData in data2:
                    residue_contact_pairs[data1].add(tuple(resData))

            for resData in residue_contact_pairs.keys():
                residue_contact_pairs[resData] = sorted(
                    residue_contact_pairs[resData])

            return residue_contact_pairs

        else:

            # get the contact atoms
            contact_atoms = self.get_contact_atoms(
                cutoff=cutoff,
                allchains=allchains,
                chain1=chain1,
                chain2=chain2,
                excludeH=excludeH,
                only_backbone_atoms=only_backbone_atoms,
                return_contact_pairs=False)

            # get the residue info
            data = dict()
            residue_contact = dict()

            for chain in contact_atoms.keys():
                data[chain] = self.get(
                    'chainID,resSeq,resName',
                    rowID=contact_atoms[chain])
                residue_contact[chain] = sorted(
                    set([tuple(resData) for resData in data[chain]]))

            return residue_contact
