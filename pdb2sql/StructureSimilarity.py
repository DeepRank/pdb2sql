import warnings
import numpy as np
from .pdb2sqlcore import pdb2sql
from .interface import interface
from .superpose import get_trans_vect, get_rotation_matrix, superpose_selection

from . import transform
import os
import pickle


class StructureSimilarity(object):

    def __init__(self, decoy, ref, verbose=False, enforce_residue_matching=True):
        """Compute structure similarity between two structures.

        This class allows to compute the i-RMSD, L-RMSD, Fnat and DockQ
        score of a given conformation.
        This can be a replacement for ProFIT.
        Note that the calculation of the zones are done by the class
        itself and does not require any extra input.

        Note:
            1. The decoy and pdb must have consistent residue numbering.
            2. The lzone files here are different with those from ProFit.
                lzone: here need only zone residues for fitting, no need
                of residue for rms calculation. RMS residues are
                automatically assumed as the other chain,
                Be careful with ProFit zone files that contain RZONE/RATOMS.
            3. Missing residues/atoms will be ignored.

        Args:
            decoy : pdb file or sql database of the decoy conformation
            ref : pdb file or sql database of the reference conformation
            verbose (bool) : verbosity option

        Examples:
            >>> from pdb2sql import StructureSimilarity
            >>> decoy = '1AK4_5w.pdb'
            >>> ref = '1AK4.pdb'
            >>> sim = StructureSimilarity(decoy,ref)
            >>> irmsd_fast = sim.compute_irmsd_fast(method='svd',
            ...     izone='1AK4.izone')
            >>> irmsd = sim.compute_irmsd_pdb2sql(method='svd',
            ...     izone='1AK4.izone')
            >>> lrmsd_fast = sim.compute_lrmsd_fast(method='svd',
            ...     lzone='1AK4.lzone',check=True)
            >>> lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None,
            ...     method='svd')
            >>> Fnat = sim.compute_fnat_pdb2sql()
            >>> Fnat_fast = sim.compute_fnat_fast(
            ...     ref_pairs='1AK4.ref_pairs')
            >>> dockQ = sim.compute_DockQScore(Fnat_fast,
            ...     lrmsd_fast,irmsd_fast)
        """

        self.decoy = decoy
        self.ref = ref
        self.verbose = verbose
        self.origin = [0., 0., 0.]
        self.enforce_residue_matching = enforce_residue_matching
        
    def __repr__(self):
        return f'{self.__module__}.{self.__class__.__name__}({self.decoy}, {self.ref}, {self.verbose})'

    def check_residues(self):
        """Check if the residue numbering matches."""

        res_ref = pdb2sql(self.ref).get_residues()
        res_dec = pdb2sql(self.decoy).get_residues()

        if res_ref != res_dec:
            print('Residues are different in the reference and decoy')
            print('Residues found in %s and not in %s' %
                  (self.ref, self.decoy))
            print(set(res_ref).difference(set(res_dec)))
            print('Residues found in %s and not in %s' %
                  (self.decoy, self.ref))
            print(set(res_dec).difference(set(res_ref)))  
            
            if self.enforce_residue_matching == True:
                raise ValueError(
                    'Residue numbering not identical in ref and decoy\n Set enforce_residue_matching=False to bypass this error.')     
            else:
                warns.Warning('Residue numbering not identical in ref and decoy.')

    ##########################################################################
    #
    #   FAST ROUTINE TO COMPUTE THE L-RMSD
    #   Require the precalculation of the lzone
    #   A dedicated routine is implemented to comoute the lzone
    #   if lzone is not given in argument the routine will compute them automatically
    #
    ##########################################################################

    # compute the L-RMSD
    def compute_lrmsd_fast(self, lzone=None, method='svd', check=True, name=['C', 'CA', 'N', 'O']):
        """Fast routine to compute the L-RMSD.

        L-RMSD is computed by aligning the longest chain of the decoy to
        the one of the reference and computing the RMSD of the shortest
        chain between decoy and reference. By default, both fitting and
         rms calculation use only backbone atoms. See reference:

        DockQ: A Quality Measure for Protein-Protein Docking Models
        https://doi.org/10.1371/journal.pone.0161879

        Args:
            lzone (None, optional): name of the file containing the zone
                definition. If None the file will be calculated first.
            method (str, optional): Method to align the fragments,
                'svd' or 'quaternion'.
            check (bool, optional): Check if the sequences are aligned
                and fix it if not. Defaults to True.
            name (list, optional): atom name to include in the zone.
                                   Defaults to ['C', 'CA', 'N', 'O']

        Returns:
            float: L-RMSD value of the conformation

        See also:
            :meth:`compute_lrmsd_pdb2sql`
        """

        # create/read the lzone file
        if lzone is None:
            resData = self.compute_lzone(save_file=False)
        elif not os.path.isfile(lzone):
            resData = self.compute_lzone(
                save_file=True, filename=lzone)
        else:
            resData = self.read_zone(lzone)

        if check:

            # Note:
            # 1. get_data_zone_backbone returns in_zone and not_in_zone
            #  here the in_zone defines the zone for fitting,
            #  and not_in_zone defines the zone for rms calculation.

            self.check_residues()

            data_decoy_long, data_decoy_short = self.get_data_zone_backbone(
                self.decoy, resData, return_not_in_zone=True, name=name)

            data_ref_long, data_ref_short = self.get_data_zone_backbone(
                self.ref, resData, return_not_in_zone=True, name=name)

            atom_long = data_ref_long.intersection(data_decoy_long)
            xyz_decoy_long = self._get_xyz(self.decoy, atom_long)
            xyz_ref_long = self._get_xyz(self.ref, atom_long)

            if data_decoy_short.symmetric_difference(data_ref_short) != set():
                res = data_decoy_short.symmetric_difference(data_ref_short)
                msg = f'\n\t Atom(s) \n {res} \n are omitted in the l-rmsd calculation'
                warnings.warn(msg)

            atom_short = data_ref_short.intersection(data_decoy_short)
            xyz_decoy_short = self._get_xyz(self.decoy, atom_short)
            xyz_ref_short = self._get_xyz(self.ref, atom_short)

        # extract the xyz
        else:

            xyz_decoy_long, xyz_decoy_short = self.get_xyz_zone_backbone(
                self.decoy, resData, return_not_in_zone=True, name=name)

            xyz_ref_long, xyz_ref_short = self.get_xyz_zone_backbone(
                self.ref, resData, return_not_in_zone=True, name=name)

        xyz_decoy_short = superpose_selection(
            xyz_decoy_short, xyz_decoy_long, xyz_ref_long, method)

        # compute the RMSD
        return self.get_rmsd(xyz_decoy_short, xyz_ref_short)

    # compute the lzone file
    def compute_lzone(self, save_file=True, filename=None):
        """Compute the zone for L-RMSD calculation.

        Note:
            It only provides the zone of long chain(s) which is used for
            fitting. The zone used for calculating RMSD is defined in
            the function `compute_lrmsd_fast`.

        Args:
            save_file (bool, optional): save the zone file
            filename (str, optional): name of the file

        Returns:
            dict: definition of the zone.
        """
        sql_ref = pdb2sql(self.ref)
        chains = list(sql_ref.get_chains())
        if len(chains) != 2:
            raise ValueError(
                'exactly two chains are needed for lrmsd calculation but we found %d' % len(chains), chains)

        nA = len(sql_ref.get('x,y,z', chainID=chains[0]))
        nB = len(sql_ref.get('x,y,z', chainID=chains[1]))

        # detect which chain is the longest
        long_chain = chains[0]
        if nA < nB:
            long_chain = chains[1]

        # extract data about the residue
        data_test = [
            tuple(data) for data in sql_ref.get(
                'chainID,resSeq',
                chainID=long_chain)]

        data_test = sorted(set(data_test))

        # close the sql
        sql_ref._close()

        if save_file:
            if filename is None:
                f = open(self.ref.split('.')[0] + '.lzone', 'w')
            else:
                f = open(filename, 'w')
            for res in data_test:
                chain = res[0]
                num = res[1]
                f.write('zone %s%d-%s%d\n' % (chain, num, chain, num))
            f.close()

        resData = {}
        for res in data_test:
            chain = res[0]
            num = res[1]

            if chain not in resData.keys():
                resData[chain] = []
            resData[chain].append(num)

        return resData

    ##########################################################################
    #
    #   FAST ROUTINE TO COMPUTE THE I-RMSD
    #   Require the precalculation of the izone
    #   A dedicated routine is implemented to comoute the izone
    #   if izone is not given in argument the routine will compute them automatcally
    #
    ##########################################################################

    def compute_irmsd_fast(self, izone=None, method='svd',
                           cutoff=10, check=True):
        """Fast method to compute the i-rmsd.

        i-RMSD is computed by selecting the backbone atoms of reference
        interface that is defined as any pair of heavy atoms from two
        chains within 10Å of each other.
        Align these backbone atoms as best as possible with their
        coutner part in the decoy and compute the RMSD. See reference:

        DockQ: A Quality Measure for Protein-Protein Docking Models
        https://doi.org/10.1371/journal.pone.0161879

        Args:
            izone (None, optional): file name of the zone.
                if None the zones will be calculated automatically.
            method (str, optional): Method to align the fragments,
                'svd' or 'quaternion'.
            cutoff (float, optional): cutoff for the contact atoms
            check (bool, optional): Check if the sequences are aligned
                and fix it if not. Should be True.

        Returns:
            float: i-RMSD value of the conformation

        See also:
            :meth:`compute_irmsd_pdb2sql`
        """

        # read the izone file
        if izone is None:
            resData = self.compute_izone(cutoff, save_file=False)
        elif not os.path.isfile(izone):
            resData = self.compute_izone(
                cutoff, save_file=True, filename=izone)
        else:
            resData = self.read_zone(izone)

        if check:

            self.check_residues()

            data_decoy = self.get_data_zone_backbone(
                self.decoy, resData, return_not_in_zone=False)
            data_ref = self.get_data_zone_backbone(
                self.ref, resData, return_not_in_zone=False)

            atom_common = data_ref.intersection(data_decoy)
            xyz_contact_decoy = self._get_xyz(self.decoy, atom_common)
            xyz_contact_ref = self._get_xyz(self.ref, atom_common)

        # extract the xyz
        else:
            xyz_contact_decoy = self.get_xyz_zone_backbone(
                self.decoy, resData)
            xyz_contact_ref = self.get_xyz_zone_backbone(
                self.ref, resData)

        # superpose the fragments
        xyz_contact_decoy = superpose_selection(xyz_contact_decoy,
                                                xyz_contact_decoy,
                                                xyz_contact_ref, method)

        # return the RMSD
        return self.get_rmsd(xyz_contact_decoy, xyz_contact_ref)

    def compute_izone(self, cutoff=10.0, save_file=True, filename=None):
        """Compute the zones for i-rmsd calculationss.

        Args:
            cutoff (float, optional): cutoff for the contact atoms
            save_file (bool, optional): svae file containing the zone
            filename (str, optional): filename

        Returns:
            dict: i-zone definition
        """

        sql_ref = interface(self.ref)
        chains = list(sql_ref.get_chains())
        if len(chains) != 2:
            raise ValueError(
                'exactly two chains are needed for irmsd calculation but we found %d' % len(chains), chains)

        contact_ref = sql_ref.get_contact_atoms(
            cutoff=cutoff, extend_to_residue=True, chain1=chains[0], chain2=chains[1])

        index_contact_ref = []
        for _, v in contact_ref.items():
            index_contact_ref += v

        # get the xyz and atom identifier of the decoy contact atoms
        data_test = [tuple(data) for data in sql_ref.get(
            'chainID,resSeq',
            rowID=index_contact_ref,
            name=sql_ref.backbone_atoms)]

        data_test = sorted(set(data_test))

        # close the sql
        sql_ref._close()

        if save_file:

            if filename is None:
                f = open(self.ref.split('.')[0] + '.izone', 'w')
            else:
                f = open(filename, 'w')

            for res in data_test:
                chain = res[0]
                num = res[1]
                f.write('zone %s%d-%s%d\n' % (chain, num, chain, num))
            f.close()

        resData = {}
        for res in data_test:
            chain = res[0]
            num = res[1]

            if chain not in resData.keys():
                resData[chain] = []
            resData[chain].append(num)
        return resData

    ##########################################################################
    #
    #   ROUTINE TO COMPUTE THE fnat QUICKLY
    #
    ##########################################################################

    def compute_fnat_fast(self, cutoff=5):
        """Fast method to cmpute the FNAT of the conformation.

        Fnat is the fraction of reference interface contacts preserved
        in the interface of decoy. The interface is defined as any pair
        of heavy atoms from two chains within 5Å of each other.

        Args:
            cutoff (int, optional): cutoff for the contact atoms

        Returns:
            float: FNAT value

        Raises:
            ValueError: if the decoy file is not found

        See also:
            :meth:`compute_fnat_pdb2sql`
        """
        # compute ref residue pairs
        residue_pairs_ref = self.compute_residue_pairs_ref(
            cutoff, save_file=False)

        # create a dict of the decoy data
        data_decoy = pdb2sql.read_pdb(self.decoy)

        # read the decoy data
        residue_xyz = {}
        residue_name = {}
        for line in data_decoy:

            if line.startswith('ATOM'):

                chainID = line[21]
                if chainID == ' ':
                    chainID = line[72]

                resSeq = int(line[22:26])
                resName = line[17:20].strip()
                name = line[12:16].strip()

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                key = (chainID, resSeq, resName)

                if key not in residue_xyz.keys():
                    residue_xyz[key] = []
                    residue_name[key] = []

                # if name in ['CA','C','N','O']
                # exclude Hydrogen
                if name[0] != 'H':
                    residue_xyz[key].append([x, y, z])
                    residue_name[key].append(name)

        # loop over the residue pairs of the ref
        nCommon, nTotal = 0, 0
        for resA, resB_list in residue_pairs_ref.items():
            if resA in residue_xyz.keys():
                xyzA = residue_xyz[resA]
                for resB in resB_list:
                    if resB in residue_xyz.keys():
                        xyzB = residue_xyz[resB]
                        dist_min = np.min(np.array(
                            [np.sqrt(np.sum((np.array(p1) - np.array(p2))**2))
                             for p1 in xyzA for p2 in xyzB]))
                        if dist_min <= cutoff:
                            nCommon += 1
                    nTotal += 1
            else:
                msg = f'\t FNAT: not find residue: {resA}'
                warnings.warn(msg)

        # normalize
        return round(nCommon / nTotal, 6)

    # compute the residue pair of the reference
    def compute_residue_pairs_ref(
            self,
            cutoff=5.0,
            save_file=True,
            filename=None):
        """Compute the residue pair on the reference conformation.

        Args:
            cutoff (float, optional): cutoff for the contact atoms
            save_file (bool, optional): save the file containing the
                residue pairs
            filename (None, optional): filename

        Returns:
            dict: defintition of the residue pairs
        """

        sql_ref = interface(self.ref)
        chains = list(sql_ref.get_chains())
        if len(chains) != 2:
            raise ValueError(
                'exactly two chains are needed for fnat calculation but we found %d' % len(chains), chains)
        residue_pairs_ref = sql_ref.get_contact_residues(
            cutoff=cutoff, return_contact_pairs=True, excludeH=True,
            chain1=chains[0], chain2=chains[1])
        sql_ref._close()

        if save_file:
            if filename is None:
                f = open(
                    self.ref.split('.')[0] +
                    'residue_contact_pairs.pckl',
                    'wb')
            else:
                f = open(filename, 'wb')

            # save as pickle
            pickle.dump(residue_pairs_ref, f)
            f.close()

        return residue_pairs_ref

    ##########################################################################
    #
    #   ROUTINE TO COMPUTE THE L-RMSD USING PDB2SQL
    #   DOES NOT REQUIRE THE PRECALCULATION OF ANYTHONG
    #   CAN OUTPUT THE SUPERIMPOSED STRUCTURES
    #   MUCH SLOWER THAN THE FAST ROUTINES BUT EASIER TO USE
    #
    ##########################################################################

    # compute the L-RMSD
    def compute_lrmsd_pdb2sql(self, exportpath=None, method='svd', **kwargs):
        """Slow routine to compute the L-RMSD.

        L-RMSD is computed by aligning the longest chain of the decoy to
        the one of the reference and computing the RMSD of the shortest
        chain between decoy and reference. Both fitting and rms calculation
        use only backbone atoms. See reference:

        DockQ: A Quality Measure for Protein-Protein Docking Models
        https://doi.org/10.1371/journal.pone.0161879

        Args:
            exportpath (str, optional): file name where the aligned pdbs
                are exported.
            method (str, optional): Method to align the fragments,
            'svd' or 'quaternion'.

        Kwargs: selection keywords used in the pdb2sql.get() method :
                'rowID', 'serial', 'name', 'altLoc',
                'resName', 'resSeq', 'iCode',
                'x', 'y', 'z', 'occ', 'temp', 'element', 'model'


        Returns:
            float: L-RMSD value of the conformation

        See also:
            :meth:`compute_lrmsd_fast`
        """
        backbone = ['CA', 'C', 'N', 'O']
        if 'name' not in kwargs:
            kwargs['name'] = backbone

        if 'chainID' in kwargs:
            raise ValueError(
                'do not specify chainID in compute_lrmsd_pdb2sql')

        # create the sql
        sql_decoy = pdb2sql(self.decoy, sqlfile='decoy.db')
        sql_ref = pdb2sql(self.ref, sqlfile='ref.db')

        # get the chains
        chains_decoy = sql_decoy.get_chains()
        chains_ref = sql_ref.get_chains()

        if chains_decoy != chains_ref:
            raise ValueError(
                'Chains are different in decoy and reference structure')

        chain1 = chains_decoy[0]
        chain2 = chains_decoy[1]

        # extract the pos of chains A
        xyz_decoy_A = np.array(
            sql_decoy.get('x,y,z', chainID=chain1, **kwargs))
        xyz_ref_A = np.array(sql_ref.get(
            'x,y,z', chainID=chain1, **kwargs))

        # extract the pos of chains B
        xyz_decoy_B = np.array(
            sql_decoy.get('x,y,z', chainID=chain2, **kwargs))
        xyz_ref_B = np.array(sql_ref.get(
            'x,y,z', chainID=chain2, **kwargs))

        # check the lengthes
        if len(xyz_decoy_A) != len(xyz_ref_A):
            xyz_decoy_A, xyz_ref_A = self.get_identical_atoms(
                sql_decoy, sql_ref, chain1, **kwargs)

        if len(xyz_decoy_B) != len(xyz_ref_B):
            xyz_decoy_B, xyz_ref_B = self.get_identical_atoms(
                sql_decoy, sql_ref, **kwargs)

        # detect which chain is the longest
        nA, nB = len(xyz_decoy_A), len(xyz_decoy_B)
        if nA > nB:
            xyz_decoy_long = xyz_decoy_A
            xyz_ref_long = xyz_ref_A

            xyz_decoy_short = xyz_decoy_B
            xyz_ref_short = xyz_ref_B

        else:
            xyz_decoy_long = xyz_decoy_B
            xyz_ref_long = xyz_ref_B

            xyz_decoy_short = xyz_decoy_A
            xyz_ref_short = xyz_ref_A

        # get the translation so that both A chains are centered
        tr_decoy = get_trans_vect(xyz_decoy_long)
        tr_ref = get_trans_vect(xyz_ref_long)

        # translate everything for 1
        xyz_decoy_short += tr_decoy
        xyz_decoy_long += tr_decoy

        # translate everuthing for 2
        xyz_ref_short += tr_ref
        xyz_ref_long += tr_ref

        # get the ideal rotation matrix
        # to superimpose the A chains
        U = get_rotation_matrix(
            xyz_decoy_long, xyz_ref_long, method=method)

        # rotate the entire fragment
        xyz_decoy_short = transform.rotate(
            xyz_decoy_short, U, center=self.origin)

        # compute the RMSD
        lrmsd = self.get_rmsd(xyz_decoy_short, xyz_ref_short)

        # export the pdb for verifiactions
        if exportpath is not None:

            # extract the pos of the dimer
            xyz_decoy = np.array(sql_decoy.get('x,y,z'))
            xyz_ref = np.array(sql_ref.get('x,y,z'))

            # translate
            xyz_ref += tr_ref
            xyz_decoy += tr_decoy

            # rotate decoy
            xyz_decoy = transform.rotate(
                xyz_decoy, U, center=self.origin)

            # update the sql database
            sql_decoy.update_column('x', xyz_decoy[:, 0])
            sql_decoy.update_column('y', xyz_decoy[:, 1])
            sql_decoy.update_column('z', xyz_decoy[:, 2])

            sql_ref.update_column('x', xyz_ref[:, 0])
            sql_ref.update_column('y', xyz_ref[:, 1])
            sql_ref.update_column('z', xyz_ref[:, 2])

            # export
            sql_decoy.exportpdb(exportpath + '/lrmsd_decoy.pdb')
            sql_ref.exportpdb(exportpath + '/lrmsd_ref.pdb')

        # close the db
        sql_decoy._close()
        sql_ref._close()

        return lrmsd

    # RETURN THE ATOMS THAT ARE SHARED BY THE TWO DB
    # FOR A GIVEN CHAINID
    @staticmethod
    def get_identical_atoms(db1, db2, chain, **kwargs):
        """Return that atoms shared by both databse for a specific chain.

        Args:
            db1 (TYPE): pdb2sql database of the first conformation
            db2 (TYPE): pdb2sql database of the 2nd conformation
            chain (str): chain name

        Kwargs: selection keywords used in the pdb2sql.get() method :
                'rowID', 'serial', 'name', 'altLoc',
                'resName', 'chainID', 'resSeq', 'iCode',
                'x', 'y', 'z', 'occ', 'temp', 'element', 'model'

        Returns:
            list, list: list of xyz for both database
        """

        # get data
        data1 = db1.get('chainID,resSeq,name',
                        chainID=chain, **kwargs)
        data2 = db2.get('chainID,resSeq,name',
                        chainID=chain, **kwargs)

        # tuplify
        data1 = [tuple(d1) for d1 in data1]
        data2 = [tuple(d2) for d2 in data2]

        # get the intersection
        shared_data = list(set(data1).intersection(data2))

        # get the xyz
        xyz1, xyz2 = [], []
        for data in shared_data:
            query = 'SELECT x,y,z from ATOM WHERE chainID=? AND resSeq=? and name=?'
            xyz1.append(list(list(db1.c.execute(query, data))[0]))
            xyz2.append(list(list(db2.c.execute(query, data))[0]))

        return xyz1, xyz2

    ##########################################################################
    #
    #   ROUTINE TO COMPUTE THE I-RMSD USING PDB2SQL
    #   DOES NOT REQUIRE THE PRECALCULATION OF ANYTHiNG
    #   BUT CAN READ AN IZONE FILE AS WELL
    #   CAN OUTPUT THE SUPERIMPOSED STRUCTURES
    #   MUCH SLOWER THAN THE FAST ROUTINES BUT EASIER TO USE
    #
    ##########################################################################

    def compute_irmsd_pdb2sql(
            self,
            cutoff=10,
            method='svd',
            izone=None,
            exportpath=None):
        """Slow method to compute the i-rmsd.

        i-RMSD is computed by selecting the backbone atoms of reference
        interface that is defined as any pair of heavy atoms from two
        chains within 10Å of each other.
        Align these backbone atoms as best as possible with their
        coutner part in the decoy and compute the RMSD. See reference:

        DockQ: A Quality Measure for Protein-Protein Docking Models
        https://doi.org/10.1371/journal.pone.0161879

        Args:
            izone (None, optional): file name of the zone.
                if None the zones will be calculated first.
            method (str, optional): Method to align the fragments,
                'svd' or 'quaternion'.
            cutoff (float, optional): cutoff for the contact atoms
            exportpath (str, optional): file name where the aligned pdbs
                are exported.

        Returns:
            float: i-RMSD value of the conformation

        See also:
            :meth:`compute_irmsd_fast`
        """

        # create thes sql
        sql_decoy = interface(self.decoy)
        sql_ref = interface(self.ref)

        # get the chains
        chains_decoy = sql_decoy.get_chains()
        chains_ref = sql_ref.get_chains()

        if chains_decoy != chains_ref:
            raise ValueError(
                'Chains are different in decoy and reference structure')

        # get the contact atoms
        if izone is None:

            contact_ref = sql_ref.get_contact_atoms(
                cutoff=cutoff,
                extend_to_residue=True,
                chain1=chains_ref[0],
                chain2=chains_ref[1])

            index_contact_ref = []
            for v in contact_ref.values():
                index_contact_ref += v
            index_contact_ref = sql_ref.get(
                'rowID', rowID=index_contact_ref, name=sql_ref.backbone_atoms)
        else:
            index_contact_ref = self.get_izone_rowID(
                sql_ref, izone, return_only_backbone_atoms=True)

        # get the xyz and atom identifier of the decoy contact atoms
        xyz_contact_ref = sql_ref.get(
            'x,y,z', rowID=index_contact_ref)
        data_contact_ref = sql_ref.get(
            'chainID,resSeq,resName,name',
            rowID=index_contact_ref)

        # get the xyz and atom indeitifier of the reference
        xyz_decoy = sql_decoy.get('x,y,z')
        data_decoy = sql_decoy.get('chainID,resSeq,resName,name')

        # loop through the ref label
        # check if the atom is in the decoy
        # if yes -> add xyz to xyz_contact_decoy
        # if no  -> remove the corresponding to xyz_contact_ref
        xyz_contact_decoy = []
        index_contact_decoy = []
        clean_ref = False
        for iat, atom in enumerate(data_contact_ref):

            try:
                index = data_decoy.index(atom)
                index_contact_decoy.append(index)
                xyz_contact_decoy.append(xyz_decoy[index])
            except Exception:
                xyz_contact_ref[iat] = None
                index_contact_ref[iat] = None
                clean_ref = True

        # clean the xyz
        if clean_ref:
            xyz_contact_ref = [
                xyz for xyz in xyz_contact_ref if xyz is not None]
            index_contact_ref = [
                ind for ind in index_contact_ref if ind is not None]

        # check that we still have atoms in both chains
        chain_decoy = list(
            set(sql_decoy.get('chainID', rowID=index_contact_decoy)))
        chain_ref = list(
            set(sql_ref.get('chainID', rowID=index_contact_ref)))

        if len(chain_decoy) < 1 or len(chain_ref) < 1:
            raise ValueError(
                'Error in i-rmsd: only one chain represented in one chain')

        # get the translation so that both A chains are centered
        tr_decoy = get_trans_vect(xyz_contact_decoy)
        tr_ref = get_trans_vect(xyz_contact_ref)

        # translate everything
        xyz_contact_decoy += tr_decoy
        xyz_contact_ref += tr_ref

        # get the ideql rotation matrix
        # to superimpose the A chains
        rot_mat = get_rotation_matrix(
            xyz_contact_decoy,
            xyz_contact_ref,
            method=method)

        # rotate the entire fragment
        xyz_contact_decoy = transform.rotate(
            xyz_contact_decoy, rot_mat, center=self.origin)

        # compute the RMSD
        irmsd = self.get_rmsd(xyz_contact_decoy, xyz_contact_ref)

        # export the pdb for verifiactions
        if exportpath is not None:

            # update the sql database
            sql_decoy.update_xyz(
                xyz_contact_decoy, rowID=index_contact_decoy)
            sql_ref.update_xyz(
                xyz_contact_ref, rowID=index_contact_ref)

            sql_decoy.exportpdb(
                exportpath + '/irmsd_decoy.pdb',
                rowID=index_contact_decoy)
            sql_ref.exportpdb(
                exportpath + '/irmsd_ref.pdb',
                rowID=index_contact_ref)

        # close the db
        sql_decoy._close()
        sql_ref._close()

        return irmsd

    # get the rowID of all the atoms
    def get_izone_rowID(self, sql, izone, return_only_backbone_atoms=True):
        """Compute the index of the izone atoms.

        Args:
            sql (pdb2sql): database of the conformation
            izone (str): filename to store the zone
            return_only_backbone_atoms (bool, optional): Returns only
                the backbone atoms

        Returns:
            lis(int): index of the atoms in the zone

        Raises:
            FileNotFoundError: if the izone file is not found
        """
        # read the file
        if not os.path.isfile(izone):
            raise FileNotFoundError('i-zone file not found', izone)

        with open(izone, 'r') as f:
            data = f.readlines()

        # get the data out of it
        resData = {}
        for line in data:

            res = line.split()[1].split('-')[0]
            chainID, resSeq = res[0], int(res[1:])

            if chainID not in resData.keys():
                resData[chainID] = []

            resData[chainID].append(resSeq)

        # get the rowID
        index_contact = []

        for chainID, resSeq in resData.items():
            if return_only_backbone_atoms:
                index_contact += sql.get('rowID',
                                         chainID=chainID,
                                         resSeq=resSeq,
                                         name=['C',
                                               'CA',
                                               'N',
                                               'O'])
            else:
                index_contact += sql.get('rowID',
                                         chainID=chainID, resSeq=resSeq)

        return index_contact

    ##########################################################################
    #
    #   ROUTINE TO COMPUTE THE fnat USING PDB2SQL
    #
    ##########################################################################

    def compute_fnat_pdb2sql(self, cutoff=5.0):
        """Slow method to compute the FNAT of the conformation.

        Fnat is the fraction of reference interface contacts preserved
        in the interface of decoy. The interface is defined as any pair
        of heavy atoms from two chains within 5Å of each other.

        Args:
            cutoff (int, optional): cutoff for the contact atoms

        Returns:
            float: FNAT value

        See also:
            :meth:`compute_fnat_fast`
        """

        # create the sql
        sql_decoy = interface(self.decoy, fix_chainID=True)
        sql_ref = interface(self.ref, fix_chainID=True)
        chains = list(sql_ref.get_chains())
        if len(chains) != 2:
            raise ValueError(
                'exactly two chains are needed for irmsd calculation but we found %d' % len(chains), chains)

        # get the contact atoms
        residue_pairs_decoy = sql_decoy.get_contact_residues(
            cutoff=cutoff, return_contact_pairs=True, excludeH=True,
            chain1=chains[0], chain2=chains[1])
        residue_pairs_ref = sql_ref.get_contact_residues(
            cutoff=cutoff, return_contact_pairs=True, excludeH=True,
            chain1=chains[0], chain2=chains[1])

        # form the pair data
        data_pair_decoy = []
        for resA, resB_list in residue_pairs_decoy.items():
            data_pair_decoy += [(resA, resB) for resB in resB_list]

        # form the pair data
        data_pair_ref = []
        for resA, resB_list in residue_pairs_ref.items():
            data_pair_ref += [(resA, resB) for resB in resB_list]

        # find the umber of residue that ref and decoys hace in common
        nCommon = len(
            set(data_pair_ref).intersection(data_pair_decoy))

        # normalize
        fnat = nCommon / len(data_pair_ref)

        sql_decoy._close()
        sql_ref._close()

        return round(fnat, 6)

    ##########################################################################
    #
    #   HELPER ROUTINES TO HANDLE THE ZONE FILES
    #
    ##########################################################################
    @staticmethod
    def get_xyz_zone_backbone(pdb_file, resData, return_not_in_zone=False, name=['C', 'CA', 'N', 'O']):
        """Get the xyz of zone backbone atoms.

        Args:
            pdb_file (str): filename containing the pdb of the molecule
            resData (dict): information about the zone residues
            return_not_in_zone (bool, optional): Do we return the
                backbone atoms not in the zone and the chains used
                in the zone.

        Returns:
            list(float): XYZ of of backbone atoms in the zone.
        """

        # read the ref file
        data = pdb2sql.read_pdb(pdb_file)

        # get the xyz of the
        xyz_in_zone = []
        xyz_not_in_zone = []

        for line in data:
            if line.startswith('ATOM'):
                chainID = line[21]
                if chainID == ' ':
                    chainID = line[72]

                resSeq = int(line[22:26])
                atname = line[12:16].strip()

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                if atname in name:
                    if chainID in resData.keys():
                        if resSeq in resData[chainID]:
                            xyz_in_zone.append([x, y, z])
                    else:
                        xyz_not_in_zone.append([x, y, z])

        if return_not_in_zone:
            return xyz_in_zone, xyz_not_in_zone

        else:
            return xyz_in_zone

    @staticmethod
    def get_data_zone_backbone(pdb_file, resData, return_not_in_zone=False, name=['C', 'CA', 'N', 'O']):
        """Get the data (chainID, resSeq, name) of backbone atoms in the zone.

        Args:
            pdb_file (str): filename containing the pdb of the molecule
            resData (dict): information about the zone residues
            return_not_in_zone (bool, optional): Do we return the atoms
                not in the zone and the chains used in the zone

        Returns:
            set(float): data of the backbone atoms in the zone
        """
        # read the ref file
        data = pdb2sql.read_pdb(pdb_file)

        # get the xyz of the
        data_in_zone = []
        data_not_in_zone = []

        for line in data:

            if line.startswith('ATOM'):

                chainID = line[21]
                if chainID == ' ':
                    chainID = line[72]

                resSeq = int(line[22:26])
                atname = line[12:16].strip()

                if atname in name:
                    if chainID in resData.keys():
                        if resSeq in resData[chainID]:
                            data_in_zone.append((chainID, resSeq, atname))
                    else:
                        data_not_in_zone.append((chainID, resSeq, atname))

        if return_not_in_zone:
            return set(data_in_zone), set(data_not_in_zone)

        else:
            return set(data_in_zone)

    @staticmethod
    def read_zone(zone_file):
        """Read the zone file.

        Args:
            zone_file (str): name of the file

        Returns:
            dict: Info about the residues in the zone

        Raises:
            FileNotFoundError: if the zone file is not found
        """
        # read the izone file
        if not os.path.isfile(zone_file):
            raise FileNotFoundError('zone file not found', zone_file)

        with open(zone_file, 'r') as f:
            data = f.readlines()

        # get the data out of it
        resData = {}
        for line in data:
            # line = zone A4-A4   for positive resNum
            # or line = zone A-4-A-4 for negative resNum
            # that happens for example in 2OUL

            # split the line
            res = line.split()[1].split('-')

            # if the resnum was positive
            # we have e.g res = [A4,A4]
            if len(res) == 2:
                res = res[0]
                chainID, resSeq = res[0], int(res[1:])

            # if the resnum was negative was negtive
            # we have e.g res = [A,4,A,4]
            elif len(res) == 4:
                chainID, resSeq = res[0], -int(res[1])

            if chainID not in resData.keys():
                resData[chainID] = []

            resData[chainID].append(resSeq)

        return resData

    @staticmethod
    def _get_xyz(pdb_file, index):
        """Get xyz using (chainID, resSeq, name) index.

        Args:
            pdb_file(file): pdb file or data
            index(set): set of index represeneted with (chainID, resSeq, name)

        Returns:
            list: list of xyz
        """
        data = pdb2sql.read_pdb(pdb_file)
        xyz = []

        for line in data:
            if line.startswith('ATOM'):
                chainID = line[21]
                if chainID == ' ':
                    chainID = line[72]

                resSeq = int(line[22:26])
                name = line[12:16].strip()

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                if (chainID, resSeq, name) in index:
                    xyz.append([x, y, z])

        return xyz

    ##########################################################################
    #
    #   CAPRI categories and DockQ score
    #
    ##########################################################################
    @staticmethod
    def compute_CapriClass(fnat, lrmsd, irmsd, system='protein-protein'):
        """Compute CAPRI ranking classes.

        Note:
            Criteria of CAPRI classes:
            https://doi.org/10.1371/journal.pone.0161879
            https://doi.org/10.1002/prot.21804
            The protocol for classifying predicted model into the four CAPRI
            categories should start with those defining incorrect predictions.

        Args:
            fnat(float): fnat
            lrmsd(float): ligand rmsd
            irmsd(float ): interface rmsd
            system (str): the type of complex system.
                Defaults to 'protein-protein'.

        Returns:
            str: CAPRI rank class, i.e. high, medium, acceptable or incorrect.
        """

        if system == 'protein-protein':
            if fnat < 0.1 or (lrmsd > 10.0 and irmsd > 4.0):
                label = 'incorrect'
            elif 0.1 <= fnat < 0.3 and (lrmsd <= 10.0 or irmsd <= 4.0) or \
                (fnat >= 0.3 and lrmsd > 5.0 and irmsd > 2.0):
                label = 'acceptable'
            elif 0.3 <= fnat < 0.5 and (lrmsd <= 5.0 or irmsd <= 2.0) or \
                (fnat >= 0.5 and lrmsd > 1.0 and irmsd > 1.0):
                label = 'medium'
            elif fnat >= 0.5 and (lrmsd <= 1.0 or irmsd <= 1.0):
                label = 'high'
        else:
            warnings.warn(
                f'Invalid complex type {system} for CAPRI class calculation')

        return label

    # compute the DockQ score from the different elements
    @staticmethod
    def compute_DockQScore(fnat, lrmsd, irmsd, d1=8.5, d2=1.5):
        """Compute the DockQ Score.

        Args:
            Fnat (float): Fnat value
            lrmsd (float): lrmsd value
            irmsd (float): irmsd value
            d1 (float, optional): first coefficient for the DockQ
                calculations
            d2 (float, optional): second coefficient for the DockQ
                calculations

        Returns:
            float: dockQ value
        """

        def scale_rms(rms, d):
            return(1. / (1 + (rms / d)**2))

        dockq = 1. / 3 * \
            (fnat + scale_rms(lrmsd, d1) + scale_rms(irmsd, d2))
        return round(dockq, 6)

    ##########################################################################
    #
    #   clahses
    #
    ##########################################################################

    @staticmethod
    def compute_clashes(pdb, chain1='A', chain2='B'):
        """Compute number of atomic clashes.

        Note:
            Clashes were defined as contacts between nonhydrogen atoms
            separated by <3.0Å. Structural models where number of clashes
            was 2 SD away from the average are excluded for assessment in
            CAPRI. see ref:  https://doi.org/10.1002/prot.10393

        Args:
            pdb(file): pdb file or data
            chain1 (str): first chain ID. Defaults to 'A'.
            chain2 (str): second chain ID. Defaults to 'B'.

        Returns:
            int: number of atomic clashes.
        """
        db = interface(pdb)
        atom_contact_pairs = db.get_contact_atoms(
            cutoff=3.0, excludeH=True,
            return_contact_pairs = True,
            chain1=chain1, chain2=chain2)
        db._close()
        nclash = 0
        for v in atom_contact_pairs.values():
            nclash += len(v)
        return nclash

    ##########################################################################
    #
    #   ROUTINES TO ACTUALY SUPERPOSE THE MOLECULES
    #
    ##########################################################################

    # compute the RMSD of two sets of points
    @staticmethod
    def get_rmsd(P, Q):
        """compute the RMSD.

        Args:
            P (np.array(nx3)): position of the points in the first
                molecule
            Q (np.array(nx3)): position of the points in the second
                molecule

        Returns:
            float: RMSD value
        """
        n = len(P)
        return round(np.sqrt(1. / n * np.sum((P - Q)**2)), 3)
