import sqlite3
import subprocess as sp
import os
import numpy as np
from time import time


class pdb2sql_base(object):

    def __init__(
            self,
            pdbfile,
            sqlfile=None,
            fix_chainID=False,
            verbose=False,
            no_extra=True):
        '''Base class for the definition of sql database.

        Args:
            pdbbfile : name of the pdbfile
            sqlfile : name of the sql file (if None the db is stored in memeory)
            fix_chainID : bool to rename chain ID from A, B, C, ....
            verbose : bool verbose
            no_extra : bool don't consider the 'temp' and 'model' column
        '''

        self.pdbfile = pdbfile
        self.sqlfile = sqlfile
        self.is_valid = True
        self.verbose = verbose
        self.no_extra = no_extra

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

        # get the data
        data = self.get('*', **kwargs)

        # write each line
        # the PDB format is pretty strict
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        if append:
            f = open(fname, 'a')
        else:
            f = open(fname, 'w')

        for d in data:
            line = 'ATOM  '
            line += '{:>5}'.format(d[0])    # serial
            line += ' '
            line += '{:^4}'.format(d[1])    # name
            line += '{:>1}'.format(d[2])    # altLoc
            line += '{:>3}'.format(d[3])  # resname
            line += ' '
            line += '{:>1}'.format(d[4])    # chainID
            line += '{:>4}'.format(d[5])    # resSeq
            line += '{:>1}'.format(d[6])    # iCODE
            line += '   '
            line += '{: 8.3f}'.format(d[7])  # x
            line += '{: 8.3f}'.format(d[8])  # y
            line += '{: 8.3f}'.format(d[9])  # z
            if not self.no_extra:
                line += '{: 6.2f}'.format(d[10])    # occ
                line += '{: 6.2f}'.format(d[11])    # temp
                line += '{:>2}'.format(d[12])       # element
            line += '\n'

            f.write(line)

        # close
        f.close()

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
