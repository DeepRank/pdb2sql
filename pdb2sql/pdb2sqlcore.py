import sqlite3
import subprocess as sp
import os
import numpy as np
from time import time

from .pdb2sql_base import pdb2sql_base


class pdb2sql(pdb2sql_base):

    def __init__(
            self,
            pdbfile,
            sqlfile=None,
            fix_chainID=False,
            verbose=False):

        super().__init__(pdbfile, sqlfile, fix_chainID, verbose)

        # create the database
        self._create_sql()

        # fix the chain ID
        if fix_chainID:
            self._fix_chainID()

        # hard limit for the number of SQL varaibles
        self.SQLITE_LIMIT_VARIABLE_NUMBER = 999
        self.max_sql_values = 950

    def _create_sql(self):
        '''Create a sql database containg a model PDB.'''
        pdbfile = self.pdbfile
        sqlfile = self.sqlfile

        if self.verbose:
            print('-- Create SQLite3 database')

         # name of the table
        table = 'ATOM'

        # size of the things
        ncol = len(self.col)
        ndel = len(self.delimiter)

        # open the data base
        # if we do not specify a db name
        # the db is only in RAM
        if self.sqlfile is None:
            self.conn = sqlite3.connect(':memory:')

        # or we create a new db file
        else:
            if os.path.isfile(sqlfile):
                sp.call('rm %s' % sqlfile, shell=True)
            self.conn = sqlite3.connect(sqlfile)
        self.c = self.conn.cursor()

        # intialize the header/placeholder
        header, qm = '', ''
        for ic, (colname, coltype) in enumerate(self.col.items()):
            header += '{cn} {ct}'.format(cn=colname, ct=coltype)
            qm += '?'
            if ic < ncol - 1:
                header += ', '
                qm += ','

        # create the table
        query = 'CREATE TABLE ATOM ({hd})'.format(hd=header)
        self.c.execute(query)

        # get pdb data
        pdbdata = pdb2sql.read_pdb(pdbfile)

        self.nModel = 0
        _check_format_ = True
        data_atom = []

        for line in pdbdata:
            line = str(line)

            if line.startswith('ATOM'):
                line = line.split('\n')[0]

            elif line.startswith('ENDMDL'):
                self.nModel += 1
                continue

            else:
                continue

            # old format chain ID fix
            if _check_format_:
                del_copy = self.delimiter.copy()
                if line[del_copy['chainID'][0]] == ' ':
                    del_copy['chainID'] = [72, 73]
                if line[del_copy['chainID'][0]] == ' ':
                    raise ValueError('chainID not found sorry')
                _check_format_ = False

            # browse all attribute of each atom
            at = ()
            for ik, (colname, coltype) in enumerate(self.col.items()):

                # get the piece of data
                if colname in del_copy.keys():
                    data = line[del_copy[colname][0]:del_copy[colname][1]].strip()

                    # convert it if necessary
                    if coltype == 'INT':
                        data = int(data)
                    elif coltype == 'REAL':
                        data = float(data)

                    # get element if it does not exist
                    if colname == "element" and not data:
                        data = pdb2sql._get_element(line)

                    # append keep the comma !!
                    # we need proper tuple
                    at += (data,)

            # append the model number
            at += (self.nModel,)

            # append
            data_atom.append(at)

        # push in the database
        self.c.executemany(
            'INSERT INTO ATOM VALUES ({qm})'.format(
                qm=qm), data_atom)

    @staticmethod
    def read_pdb(pdbfile):
        # read the pdb file a pure python way
        # RMK we go through the data twice here. Once to read the ATOM line and once to parse the data ...
        # we could do better than that. But the most time consuming step seems
        # to be the CREATE TABLE query
        if isinstance(pdbfile, str):
            if os.path.isfile(pdbfile):
                with open(pdbfile, 'r') as fi:
                    pdbdata = fi.readlines()
            else:
                raise FileNotFoundError(f'PDB file {pdbfile} not found')
        elif isinstance(pdbfile, list):
            if isinstance(pdbfile[0], str):
                pdbdata = pdbfile
            elif isinstance(pdbfile[0], bytes):
                pdbdata = [line.decode() for line in pdbfile]
            else:
                raise ValueError(f'Non-valid pdb input {pdbfile}')
        elif isinstance(pdbfile, np.ndarray):
            pdbfile = pdbfile.tolist()
            if isinstance(pdbfile[0], str):
                pdbdata = pdbfile
            elif isinstance(pdbfile[0], bytes):
                pdbdata = [line.decode() for line in pdbfile]
            else:
                raise ValueError(f'Non-valid pdb input {pdbfile}')
        else:
            raise ValueError(f'Non-valid pdb input: {pdbfile}')

        return pdbdata

    @staticmethod
    def _get_element(pdb_line):
        """Get element type from the atom type of a pdb line

        Notes:
            Atom type occupies 13-16th columns of a PDB line.
            http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            Four situations exist:
                13 14 15 16
                   C  A      The element is C
                C  A         The element is Ca
                1  H  G      The element is H
                H  E  2  1   The element is H

        Args:
            pdb_line(str): one PDB ATOM line

        Returns:
            [str]: element name
        """
        first_char = pdb_line[12].strip()
        last_char = pdb_line[15].strip()
        if first_char:
            if first_char in "0123456789":
                elem = pdb_line[13]
            elif first_char == "H" and last_char:
                elem = "H"
            else:
                elem = pdb_line[12:14]
        else:
            elem = pdb_line[13]
        return elem

    # replace the chain ID by A,B,C,D, ..... in that order
    def _fix_chainID(self):

        from string import ascii_uppercase

        # get the current names
        chainID = self.get('chainID')
        natom = len(chainID)
        chainID = sorted(set(chainID))

        if len(chainID) > 26:
            print(
                "Warning more than 26 chains have been detected. This is so far not supported")
            sys.exit()

        # declare the new names
        newID = [''] * natom

        # fill in the new names
        for ic, chain in enumerate(chainID):
            index = self.get('rowID', chain=chain)
            for ind in index:
                newID[ind] = ascii_uppercase[ic]

        # update the new name
        self.update_column('chainID', newID)

    # get the names of the columns

    def get_colnames(self):
        cd = self.conn.execute('select * from atom')
        print('Possible column names are:')
        names = list(map(lambda x: x[0], cd.description))
        print('\trowID')
        for n in names:
            print('\t' + n)

    # print the database
    def prettyprint(self):
        import pandas.io.sql as psql
        df = psql.read_sql("SELECT * FROM ATOM", self.conn)
        print(df)

    def uglyprint(self):
        ctmp = self.conn.cursor()
        ctmp.execute("SELECT * FROM ATOM")
        print(ctmp.fetchall())

    # get the properties
    def get(self, atnames, **kwargs):
        '''Exectute a simple SQL query that extracts values of attributes for certain condition.

        Args :
            attribute (str) : attribute to retreive eg : ['x','y,'z'], 'xyz', 'resSeq'
                              if None all the attributes are returned

            **kwargs : argument to select atoms eg : name = ['CA','O'], chainID = 'A', no_name = ['H']

        Returns:
            data : array containing the value of the attributes

        Example :
        >>> db.get('x,y,z',chainID='A',no_name=['H']")
        '''

        # the asked keys
        keys = kwargs.keys()

        # check if the column exists
        try:
            self.c.execute(
                "SELECT EXISTS(SELECT {an} FROM ATOM)".format(
                    an=atnames))
        except BaseException:
            print('Error column %s not found in the database' % atnames)
            self.get_colnames()
            return

        if 'model' not in kwargs.keys() and self.nModel > 0:
            model_data = []
            for iModel in range(self.nModel):
                kwargs['model'] = iModel
                model_data.append(self.get(atnames, **kwargs))
            return model_data

        # if we have 0 key we take the entire db
        if len(kwargs) == 0:
            query = 'SELECT {an} FROM ATOM'.format(an=atnames)
            data = [list(row) for row in self.c.execute(query)]

        #######################################################################
        # GENERIC QUERY
        #
        # each keys must be a valid columns
        # each valu may be a single value or an array
        # AND is assumed between different keys
        # OR is assumed for the different values of a given key
        #
        #######################################################################
        else:

            # check that all the keys exists
            for k in keys:
                if k.startswith('no_'):
                    k = k[3:]

                try:
                    self.c.execute(
                        "SELECT EXISTS(SELECT {an} FROM ATOM)".format(
                            an=k))
                except BaseException:
                    print('Error column %s not found in the database' % k)
                    self.get_colnames()
                    return

            # form the query and the tuple value
            query = 'SELECT {an} FROM ATOM WHERE '.format(an=atnames)
            conditions = []
            vals = ()

            # iterate through the kwargs
            for ik, (k, v) in enumerate(kwargs.items()):

                # deals with negative conditions
                if k.startswith('no_'):
                    k = k[3:]
                    neg = ' NOT'
                else:
                    neg = ''

                # get if we have an array or a scalar
                # and build the value tuple for the sql query
                # deal with the indexing issue if rowID is required
                if isinstance(v, list):
                    nv = len(v)

                    # if we have a large number of values
                    # we must cut that in pieces because SQL has a hard limit
                    # that is 999. The limit is here set to 950
                    # so that we can have multiple conditions with a total number
                    # of values inferior to 999
                    if nv > self.max_sql_values:

                        # cut in chunck
                        chunck_size = self.max_sql_values
                        vchunck = [v[i:i + chunck_size]
                                   for i in range(0, nv, chunck_size)]

                        data = []
                        for v in vchunck:
                            new_kwargs = kwargs.copy()
                            new_kwargs[k] = v
                            data += self.get(atnames, **new_kwargs)
                        return data

                    # otherwise we just go on
                    else:
                        if k == 'rowID':
                            vals = vals + tuple([int(iv + 1) for iv in v])
                        else:
                            vals = vals + tuple(v)
                else:
                    nv = 1
                    if k == 'rowID':
                        vals = vals + (int(v + 1),)
                    else:
                        vals = vals + (v,)

                # create the condition for that key
                conditions.append(k + neg + ' in (' + ','.join('?' * nv) + ')')

            # stitch the conditions and append to the query
            query += ' AND '.join(conditions)

            # error if vals is too long
            if len(vals) > self.SQLITE_LIMIT_VARIABLE_NUMBER:
                print('\nError : SQL Queries can only handle a total of 999 values')
                print('      : The current query has %d values' % len(vals))
                print('      : Hence it will fails.')
                print(
                    '      : You are in a rare situation where MULTIPLE conditions have')
                print('      : have a combined number of values that are too large')
                print('      : These conditions are:')
                ntot = 0
                for k, v in kwargs.items():
                    print('      : --> %10s : %d values' % (k, len(v)))
                    ntot += len(v)
                print('      : --> %10s : %d values' % ('Total', ntot))
                print('      : Try to decrease max_sql_values in pdb2sql.py\n')
                raise ValueError('Too many SQL variables')

            # query the sql database and return the answer in a list
            data = [list(row) for row in self.c.execute(query, vals)]

        # empty data
        if len(data) == 0:
            print('Warning sqldb.get returned an empty')
            return data

        # fix the python <--> sql indexes
        # if atnames == 'rowID':
        if 'rowID' in atnames:
            index = atnames.split(',').index('rowID')
            for i in range(len(data)):
                data[i][index] -= 1

        # postporcess the output of the SQl query
        # flatten it if each els is of size 1
        if len(data[0]) == 1:
            data = [d[0] for d in data]

        return data

    def update(self, attribute, values, **kwargs):
        '''Update the database.

        Args:
            attribute (str) : string of attribute names eg. ['x','y,'z'], 'xyz', 'resSeq'

            values (np.ndarray) : an array of values that corresponds
                                  to the number of attributes and atoms selected

            **kwargs : selection arguments eg : name = ['CA','O'], chainID = 'A', no_name = ['H']
        '''

        # the asked keys
        keys = kwargs.keys()

        # check if the column exists
        try:
            self.c.execute(
                "SELECT EXISTS(SELECT {an} FROM ATOM)".format(
                    an=attribute))
        except BaseException:
            print('Error column %s not found in the database' % attribute)
            self.get_colnames()
            return

        # TODO
        # handle the multi model cases
        if 'model' not in keys and self.nModel > 0:
            for iModel in range(self.nModel):
                kwargs['model'] = iModel
                self.update(attribute, values, **kwargs)
            return

        # parse the attribute
        if ',' in attribute:
            attribute = attribute.split(',')

        if not isinstance(attribute, list):
            attribute = [attribute]

        # check the size
        natt = len(attribute)
        nrow = len(values)
        ncol = len(values[0])

        if natt != ncol:
            raise ValueError(
                'Number of attribute incompatible with the number of columns in the data')

        # get the row ID of the selection
        rowID = self.get('rowID', **kwargs)
        nselect = len(rowID)

        if nselect != nrow:
            raise ValueError(
                'Number of data values incompatible with the given conditions')

        # prepare the query
        query = 'UPDATE ATOM SET '
        query = query + ', '.join(map(lambda x: x + '=?', attribute))
        query = query + ' WHERE rowID=?'

        # prepare the data
        data = []
        for i, val in enumerate(values):

            tmp_data = [v for v in val]

            # here the conversion of the indexes is a bit annoying
            tmp_data += [rowID[i] + 1]

            data.append(tmp_data)

        self.c.executemany(query, data)

    def update_column(self, colname, values, index=None):
        '''Update a single column.

        Args:
            colname (str): name of the column to update
            values (list): new values of the column
            index (None, optional): index of the column to update (default all)

        Example:
        >>> db.update_column('x',np.random.rand(10),index=list(range(10)))
        '''

        if index is None:
            data = [[v, i + 1] for i, v in enumerate(values)]
        else:
            # shouldn't that be ind+1 ?
            data = [[v, ind] for v, ind in zip(values, index)]

        query = 'UPDATE ATOM SET {cn}=? WHERE rowID=?'.format(cn=colname)
        self.c.executemany(query, data)
        # self.conn.commit()

    def add_column(self, colname, coltype='FLOAT', default=0):
        '''Add an etra column to the ATOM table.'''
        query = "ALTER TABLE ATOM ADD COLUMN '%s' %s DEFAULT %s" % (
            colname, coltype, str(default))
        self.c.execute(query)

    def commit(self):
        '''Cpmmit to the database.'''
        self.conn.commit()

    # def _create_sql_nomodel(self):
    #     '''Create the sql database.'''

    #     pdbfile = self.pdbfile
    #     sqlfile = self.sqlfile

    #     if self.verbose:
    #         print('-- Create SQLite3 database')

    #      #name of the table
    #     table = 'ATOM'

    #     # size of the things
    #     ncol = len(self.col)
    #     ndel = len(self.delimiter)

    #     # open the data base
    #     # if we do not specify a db name
    #     # the db is only in RAM
    #     if self.sqlfile is None:
    #         self.conn = sqlite3.connect(':memory:')

    #     # or we create a new db file
    #     else:
    #         if os.path.isfile(sqlfile):
    #             sp.call('rm %s' %sqlfile,shell=True)
    #         self.conn = sqlite3.connect(sqlfile)
    #     self.c = self.conn.cursor()

    #     # intialize the header/placeholder
    #     header,qm = '',''
    #     for ic,(colname,coltype) in enumerate(self.col.items()):
    #         header += '{cn} {ct}'.format(cn=colname,ct=coltype)
    #         qm += '?'
    #         if ic < ncol-1:
    #             header += ', '
    #             qm += ','

    #     # create the table
    #     query = 'CREATE TABLE ATOM ({hd})'.format(hd=header)
    #     self.c.execute(query)

    #     # read the pdb file a pure python way
    #     # RMK we go through the data twice here. Once to read the ATOM line and once to parse the data ...
    #     # we could do better than that. But the most time consuming step seems to be the CREATE TABLE query
    #     with open(pdbfile,'r') as fi:
    #         data = [line.split('\n')[0] for line in fi if line.startswith('ATOM')]

    #     # if there is no ATOM in the file
    #     if len(data)==1 and data[0]=='':
    #         print("-- Error : No ATOM in the pdb file.")
    #         self.is_valid = False
    #         return

    #     # old format chain ID fix
    #     del_copy = self.delimiter.copy()
    #     if data[0][del_copy['chainID'][0]] == ' ':
    #         del_copy['chainID'] = [72,73]

    #     # get all the data
    #     data_atom = []
    #     for iatom,atom in enumerate(data):

    #         # sometimes we still have an empty line somewhere
    #         if len(atom) == 0:
    #             continue

    #         # browse all attribute of each atom
    #         at = ()
    #         for ik,(colname,coltype) in enumerate(self.col.items()):

    #             # get the piece of data
    #             data = atom[del_copy[colname][0]:del_copy[colname][1]].strip()

    #             # convert it if necessary
    #             if coltype == 'INT':
    #                 data = int(data)
    #             elif coltype == 'REAL':
    #                 data = float(data)

    #             # append keep the comma !!
    #             # we need proper tuple
    #             at +=(data,)

    #         # append
    #         data_atom.append(at)

    #     # push in the database
    #     self.c.executemany('INSERT INTO ATOM VALUES ({qm})'.format(qm=qm),data_atom)
