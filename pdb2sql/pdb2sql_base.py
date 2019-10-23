import sqlite3
import subprocess as sp
import os
import warnings
import numpy as np
from time import time


class pdb2sql_base(object):

    def __init__(
            self,
            pdbfile,
            sqlfile=None,
            fix_chainID=False,
            verbose=False):
        '''Base class for the definition of sql database.

        Args:
            pdbfile (str, list(str/bytes), ndarray) : name of pdbfile or
                list or ndarray containing the pdb data
            sqlfile (str, optional): name of the sqlfile.
                By default it is created in memory only.
            fix_chainID (bool, optinal): check if the name of the chains
                are A,B,C, .... and fix it if not.
            verbose (bool): probably print stuff
        '''

        self.pdbfile = pdbfile
        self.sqlfile = sqlfile
        self.is_valid = True
        self.verbose = verbose

        # column names and types
        self.col = {'serial': 'INT',
                    'name': 'TEXT',
                    'altLoc': 'TEXT',
                    'resName': 'TEXT',
                    'chainID': 'TEXT',
                    'resSeq': 'INT',
                    'iCode': 'TEXT',
                    'x': 'REAL',
                    'y': 'REAL',
                    'z': 'REAL',
                    'occ': 'REAL',
                    'temp': 'REAL',
                    'element': 'TEXT',
                    'model': 'INT'}

        # delimtier of the column format
        # taken from
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        self.delimiter = {
            'serial': [6, 11],
            'name': [12, 16],
            'altLoc': [16, 17],
            'resName': [17, 20],
            'chainID': [21, 22],
            'resSeq': [22, 26],
            'iCode': [26, 27],
            'x': [30, 38],
            'y': [38, 46],
            'z': [46, 54],
            'occ': [54, 60],
            'temp': [60, 66],
            'element': [76,78]}

    ##########################################################################
    #
    #   CREATION AND PRINTING
    #
    ##########################################################################

    '''
    Main function to create the SQL data base
    '''

    def _create_sql(self):
        raise NotImplementedError()

    # get the properties
    def get(self, atnames, **kwargs):
        raise NotImplementedError()

    def get_xyz(self, **kwargs):
        '''shortcut to get the xyz.'''
        return self.get('x,y,z', **kwargs)

    def get_residues(self, **kwargs):
        '''Get the sequence of the selection.'''

        res = [tuple(x) for x in self.get('chainID,resName,resSeq', **kwargs)]
        return sorted(set(res), key=res.index)

    def get_chains(self, **kwargs):
        '''get the chain IDS.'''
        chains = self.get('chainID', **kwargs)
        return sorted(set(chains))

    def update(self, attribute, values, **kwargs):
        raise NotImplementedError()

    def update_xyz(self, xyz, **kwargs):
        '''Update the xyz coordinate.'''
        self.update('x,y,z', xyz, **kwargs)

    def update_column(self, colname, values, index=None):
        '''Update a single column.'''
        raise NotImplementedError()

    def add_column(self, colname, coltype='FLOAT', default=0):
        '''Add an etra column to the ATOM table.'''
        raise NotImplementedError()

    def exportpdb(self, fname, append=False, periodic=False, **kwargs):
        '''Export a PDB file with kwargs selection.'''

        if append:
            f = open(fname, 'a')
        else:
            f = open(fname, 'w')

        lines = self.sql2pdb(**kwargs)
        for i in lines:
            f.write(i + '\n')

        f.close()

    def sql2pdb(self, **kwargs):
        """Convert sql pdb data to PDB formatted lines

        Returns:
            list: pdb-format lines
        """
        data = self.get('*', **kwargs)
        pdb = []
        # the PDB format is pretty strict
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        for d in data:
            line = 'ATOM  '
            line += '{:>5}'.format(d[0])    # serial
            line += ' '
            line += self._format_atomname(d) # name
            line += '{:>1}'.format(d[2])    # altLoc
            line += '{:>3}'.format(d[3])    # resname
            line += ' '
            line += '{:>1}'.format(d[4])    # chainID
            line += '{:>4}'.format(d[5])    # resSeq
            line += '{:>1}'.format(d[6])    # iCODE
            line += '   '
            line += pdb2sql_base._format_xyz(d[7]) # x
            line += pdb2sql_base._format_xyz(d[8]) # y
            line += pdb2sql_base._format_xyz(d[9])  # z
            line += '{:>6.2f}'.format(d[10])    # occ
            line += '{:>6.2f}'.format(d[11])    # temp
            line += ' ' * 10
            line += '{:>2}'.format(d[12])       # element
            # line += '\n'
            pdb.append(line)

        return pdb

    def _format_atomname(self, data):
        """Format atom name to align with PDB reqireuments:
             - alignment of one-letter atom name starts at column 14,
             - while two-letter atom name such as FE starts at column 13.

        Args:
            data(list): sql output for one pdb line

        Returns:
            str: formatted atom name
        """
        name = data[1]
        lname = len(name)
        if lname in (1, 4):
            name = '{:^4}'.format(name)
        elif lname == 2:
            if name == data[12]:  # name == element
                name = '{:<4}'.format(name)
            else:
                name = '{:^4}'.format(name)
        else:
            if name[0] in '0123456789':
                name = '{:<4}'.format(name)
            else:
                name = '{:>4}'.format(name)
        return name

    @staticmethod
    def _format_xyz(i):
        """Format PDB coordinations x,y or z value.

        Note: PDB has a fixed 8-column space for x,y or z value.
            Thus the value should be in the range of (-1e7, 1e8).

        Args:
            i(float): PDB coordinations x, y or z.

        Raises:
            ValueError: Exceed the range of (-1e7, 1e8)

        Returns:
            str: formated x, y or z value.
        """

        if i >= 1e8 - 0.5 or i <= -1e7 + 0.5:
            raise ValueError(
                f'PDB coordination {i} exceeds the range of (-1e7, 1e8) '
                f'after rounding.')
        elif i >= 1e6 - 0.5 or i <= -1e5 + 0.5:
            i = '{:>8.0f}'.format(i)
        elif i >= 1e5 - 0.5 or i <= -1e4 + 0.5:
            i = '{:>8.1f}'.format(i)
        elif i >= 1e4 - 0.5 or i <= -1e3 + 0.5:
            i = '{:>8.2f}'.format(i)
        else:
            i = '{:>8.3f}'.format(i)

        return i

    def close(self, rmdb=True):

        if self.sqlfile is None:
            self.conn.close()

        else:

            if rmdb:
                self.conn.close()
                os.system('rm %s' % (self.sqlfile))
            else:
                self.commit()
                self.conn.close()
