import sqlite3
import warnings
import subprocess as sp
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

from .pdb2sql_base import pdb2sql_base


class pdb2sql(pdb2sql_base):

    def __init__(self, pdbfile, tablename='atom', **kwargs):
        """Create a SQL database with PDB data.

        Notes:
            Only "ATOM" data of PDB is parsed, other items e.g. HETATM
            are not ignored.

        Args:
            pdbfile(str, list, ndarray): pdb file or data

        Examples:
            >>> db = pdb2sql.pdb2sql('3CRO.pdb')
        """
        super().__init__(pdbfile, **kwargs)

        # create the database
        self._create_sql()
        self._create_table(pdbfile, tablename=tablename)

        # fix the chain ID
        if self.fix_chainID:
            self._fix_chainID()

    def __call__(self, **kwargs):
        """Returns an pdb2sql instance of the selected parts

         Args:

            kwargs: argument to select atoms, dict value must be list,
                e.g.:
                    - name = ['CA', 'O']
                    - no_name = ['CA', 'C']
                    - chainID = ['A']
                    - no_chainID = ['A']

        Returns:
            pdb2sql: an pb2sql instance

        Examples:
            >>> sqldb = pdb2sql('1AK4.pdb')
            >>> dbA = sqldb(chainID='A')
        """

        names = self._get_table_names()
        if len(names) > 1:
            warnings.warn('pdbsql is meant for single structure. \
                          To use multiple structures use many2sql. \
                          This call will only return the data of \
                          the first table : ', names[0])

        pdb_data = self.sql2pdb(tablename=names[0], **kwargs)
        new_db = pdb2sql(pdb_data, tablename=names[0])

        return new_db

    def __repr__(self):
        return f'{self.__module__}.{self.__class__.__name__} object'

    def _create_sql(self, tablename='ATOM'):
        """Create a sql database containg a model PDB."""
        sqlfile = self.sqlfile

        if self.verbose:
            print('-- Create SQLite3 database')

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

    def _create_table(self, pdbfile, tablename='ATOM'):

        # check the length of the 1st line to see if we can
        # read all data
        check_line_length = False

        # size of the things
        ncol = len(self.col)

        # clean up tablename
        for c in "!@#$%^&*()[]{};:,./<>?\|`~-=_+":
            tablename = tablename.replace(c, '_')

        # intialize the header/placeholder
        header, qm = '', ''
        for ic, (colname, coltype) in enumerate(self.col.items()):
            header += '{cn} {ct}'.format(cn=colname, ct=coltype)
            qm += '?'
            if ic < ncol - 1:
                header += ', '
                qm += ','

        # create the table
        query = 'CREATE TABLE {tablename} ({hd})'.format(tablename=tablename,
                                                         hd=header)
        self.c.execute(query)

        # get pdb data
        pdbdata = pdb2sql.read_pdb(pdbfile)

        self._nModel = 0
        data_atom = []

        for line in pdbdata:
            line = str(line)

            if line.startswith('ATOM'):
                line = line.split('\n')[0]

            elif line.startswith('ENDMDL'):
                self._nModel += 1
                continue

            else:
                continue

            if not check_line_length:
                # set the column we want to read
                self._set_col_values(line)
                check_line_length = True

            # format
            line = self._format_pdb_linelength(line)

            # browse all attribute of each atom
            at = ()
            for colname, coltype in self.col.items():

                # get the piece of data
                if colname in self.delimiter.keys():
                    data = line[self.delimiter[colname][0]:
                                self.delimiter[colname][1]].strip()

                    # check pdb format and reset values if necessary
                    # Empty chainID, occ, temp and element are not allowed
                    if not data:
                        if colname == "chainID":
                            data = pdb2sql._get_chainID(line)
                        if colname == "occ":
                            data = 1.00
                        if colname == "temp":
                            data = 10.00
                        if colname == "element":
                            data = pdb2sql._get_element(line)

                    # convert it if necessary
                    if coltype == 'INT':
                        data = int(data)
                    elif coltype == 'REAL':
                        data = float(data)

                    # append keep the comma !!
                    # we need proper tuple
                    at += (data,)

            # append the model number
            at += (self._nModel,)

            # append
            data_atom.append(at)

        # push in the database
        self.c.executemany(
            'INSERT INTO {tablename} VALUES ({qm})'.format(tablename=tablename,
                                                           qm=qm), data_atom)

    @staticmethod
    def read_pdb(pdbfile):
        """Read pdb file or data to a list.

        Args:
            pdbfile(str, list or ndarray): pdb file or data

        Raises:
            FileNotFoundError: pdb file not found
            ValueError: invalid input

        Returns:
            list: pdb content in list format

        Examples:
            >>> pdb2sql.read_pdb('3CRO.pdb')
        """
        # read the pdb file a pure python way
        # RMK we go through the data twice here. Once to read the ATOM line and once to parse the data ...
        # we could do better than that. But the most time consuming step seems
        # to be the CREATE TABLE query
        if isinstance(pdbfile, bytes):
            pdbfile = pdbfile.decode()

        if isinstance(pdbfile, str):
            if os.path.exists(pdbfile):
                if os.path.isfile(pdbfile):
                    with open(pdbfile, 'r') as fi:
                        pdbdata = fi.readlines()
                else:
                    raise FileNotFoundError(
                        f'{pdbfile} is not a file')
            else:
                # input is pdb content
                if pdbfile.count('\nATOM ') > 3:
                    pdbdata = pdbfile.split('\n')
                # invalid path
                else:
                    raise FileNotFoundError(
                        f'File not found: {pdbfile}')
        elif isinstance(pdbfile, Path):
            if not pdbfile.exists():
                raise FileNotFoundError(f'File not found: {pdbfile}')
            elif pdbfile.is_file():
                with pdbfile.open() as fi:
                    pdbdata = fi.readlines()
            else:
                raise FileNotFoundError(f'{pdbfile} is not a file')
        elif isinstance(pdbfile, list):
            if isinstance(pdbfile[0], str):
                pdbdata = pdbfile
            elif isinstance(pdbfile[0], bytes):
                pdbdata = [line.decode() for line in pdbfile]
            else:
                raise ValueError(f'Invalid pdb input: {pdbfile}')
        elif isinstance(pdbfile, np.ndarray):
            pdbfile = pdbfile.tolist()
            if isinstance(pdbfile[0], str):
                pdbdata = pdbfile
            elif isinstance(pdbfile[0], bytes):
                pdbdata = [line.decode() for line in pdbfile]
            else:
                raise ValueError(f'Invalid pdb input: {pdbfile}')
        else:
            raise ValueError(f'Invalid pdb input: {pdbfile}')

        return pdbdata

    @staticmethod
    def _format_pdb_linelength(pdb_line):
        linelen = len(pdb_line)
        if linelen < 80:
            pdb_line = pdb_line + ' ' * (80 - linelen)
        elif linelen > 80:
            pdb_line = pdb_line[:54]
        return pdb_line

    def _set_col_values(self, pdb_line):
        linelen = len(pdb_line)

        # if the line is less than 80 char
        # we read all
        if linelen < 80:
            pdb_line = pdb_line + ' ' * (80 - linelen)
            self.col = self.base_col.update(self.extra_col)
            self.delimiter = self.base_delimiter.update(self.extra_delimiter)

        # otherwise we only read the data up to the xyz col
        elif linelen > 80:
            pdb_line = pdb_line + ' ' * (80 - linelen)
            self.col = self.base_col
            self.delimiter = self.base_delimiter
            raise Warning(
                f'pdb line is longer than 80:\n{pdb_line}')

        return pdb_line

    @staticmethod
    def _get_chainID(pdb_line):
        segID_ind = [72, 76]    # segID columns in pdb
        segID = pdb_line[segID_ind[0]:segID_ind[1]].strip()
        if segID:
            warnings.warn("Missing chainID and set it with segID")
            return segID
        else:
            raise ValueError('chainID not found')

    @staticmethod
    def _get_element(pdb_line):
        """Get element type from the atom type of a pdb line.

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
            str : element name
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
        # warnings.warn(f'Element is missing and guessed using atom type for line\n {pdb_line}')
        return elem

    # replace the chain ID by A,B,C,D, ..... in that order
    def _fix_chainID(self):

        from string import ascii_uppercase

        # get the current names
        chainID = self.get('chainID')
        natom = len(chainID)
        chainID = sorted(set(chainID))

        if len(chainID) > 26:
            warnings.warn(
                "More than 26 chains have been detected, not supported so far")
            sys.exit(1)

        # declare the new names
        newID = [''] * natom

        # fill in the new names
        for ic, chain in enumerate(chainID):
            index = self.get('rowID', chainID=chain)
            for ind in index:
                newID[ind] = ascii_uppercase[ic]

        # update the new name
        self.update_column('chainID', newID)

    # get the names of the columns

    def get_colnames(self):
        """Get SQL column names.

        Returns:
            list: all available column names

        Examples:
            >>> db.get_colnames()
        """
        tablename = self._get_table_names()[0]
        cd = self.conn.execute(
            'select * from {tablename}'.format(tablename=tablename))
        names = list(map(lambda x: x[0], cd.description))
        names = ['rowID'] + names
        return names

    def print_colnames(self):
        """Print out SQL column names.

        Examples:
            >>> db.print_colnames()
        """
        tablenames = self._get_table_names()
        names = self.get_colnames()
        print('Possible column names are:')
        for n in names:
            print('\t' + n)

    def print(self, columns='*', tablename='atom', **kwargs):
        """Print out SQL ATOM table.

        Notes:
            Float number is stored in original precision.
            It will be formatted properly when output pdb with
            `exportpdb` method.

        Args:
            columns (str): columns to retreive, eg: "x,y,z".
                if "*" all the columns are returned.
                Avaiable columns:
                    serial, name, altLoc, resName, chainID, resSeq,
                    iCode, x, y, z, occ, temp, element, model

            kwargs: argument to select atoms, dict value must be list,
                e.g.:
                    - name = ['CA', 'O']
                    - no_name = ['CA', 'C']
                    - chainID = ['A']
                    - no_chainID = ['A']

        Examples:
            >>> db.print()
        """
        data = self.get(columns, tablename=tablename, **kwargs)
        arr = np.array(data)

        if len(arr.shape) == 2:
            df = pd.DataFrame(arr)
        elif len(arr.shape) == 3:
            arr = arr.reshape((-1, arr.shape[2]))
            df = pd.DataFrame(arr)

        if columns == '*':
            cd = self.conn.execute(
                'select * from {tablename}'.format(tablename=tablename))
            names = list(map(lambda x: x[0], cd.description))
            df.columns = names
        else:
            df.columns = columns.split(',')

        print(df.to_csv(sep='\t', index=False))

    # get the properties
    def get(self, columns, tablename='ATOM', **kwargs):
        """Exectute simple SQL query to extract values.

        Args:
            columns (str): columns to retreive, eg: "x,y,z".
                if "*" all the columns are returned.
                Avaiable columns:
                    serial, name, altLoc, resName, chainID, resSeq,
                    iCode, x, y, z, occ, temp, element, model

            kwargs: argument to select atoms, dict value must be list,
                e.g.:
                    - name = ['CA', 'O']
                    - no_name = ['CA', 'C']
                    - chainID = ['A']
                    - no_chainID = ['A']

        Returns:
            data: list containing the value of the attributes

        Examples:
            >>> db.get('x,y,z', chainID=['A'], no_resName=['ALA', 'TRP'])
        """

        # check arguments format
        valid_colnames = self.get_colnames()

        if not isinstance(columns, str):
            raise TypeError("argument columns must be str")

        if columns != '*':
            for i in columns.split(','):
                if i.strip() not in valid_colnames:
                    raise ValueError(
                        f'Invalid column name {i}. Possible names are\n'
                        f'{self.get_colnames()}')

        # the asked keys
        keys = kwargs.keys()

        if 'model' not in kwargs.keys() and self._nModel > 0:
            model_data = []
            for iModel in range(self._nModel):
                kwargs['model'] = iModel
                model_data.append(self.get(columns, **kwargs))
            return model_data

        # if we have 0 key we take the entire db
        if len(kwargs) == 0:
            query = 'SELECT {an} FROM {tablename}'.format(
                an=columns, tablename=tablename)
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
                        "SELECT EXISTS(SELECT {an} FROM {tablename})".format(
                            an=k, tablename=tablename))
                except BaseException:
                    raise ValueError(
                        f'Invalid column name {k}. Possible names are\n'
                        f'{self.get_colnames()}')

            # form the query and the tuple value
            query = 'SELECT {an} FROM {tablename} WHERE '.format(
                an=columns, tablename=tablename)
            conditions = []
            vals = ()

            # iterate through the kwargs
            for _, (k, v) in enumerate(kwargs.items()):

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
                            data += self.get(columns, **new_kwargs)
                        return data

                    # otherwise we just go on
                    else:
                        if k == 'rowID':
                            vals = vals + \
                                tuple([int(iv + 1) for iv in v])
                        else:
                            vals = vals + tuple(v)
                else:
                    nv = 1
                    if k == 'rowID':
                        vals = vals + (int(v + 1),)
                    else:
                        vals = vals + (v,)

                # create the condition for that key
                conditions.append(
                    k + neg + ' in (' + ','.join('?' * nv) + ')')

            # stitch the conditions and append to the query
            query += ' AND '.join(conditions)

            # error if vals is too long
            if len(vals) > self.SQLITE_LIMIT_VARIABLE_NUMBER:
                print(
                    '\nError : SQL Queries can only handle a total of 999 values')
                print('      : The current query has %d values' %
                      len(vals))
                print('      : Hence it will fails.')
                print(
                    '      : You are in a rare situation where MULTIPLE conditions have')
                print(
                    '      : have a combined number of values that are too large')
                print('      : These conditions are:')
                ntot = 0
                for k, v in kwargs.items():
                    print('      : --> %10s : %d values' %
                          (k, len(v)))
                    ntot += len(v)
                print('      : --> %10s : %d values' %
                      ('Total', ntot))
                print(
                    '      : Try to decrease max_sql_values in pdb2sql.py\n')
                raise ValueError('Too many SQL variables')

            # query the sql database and return the answer in a list
            data = [list(row) for row in self.c.execute(query, vals)]

        # empty data
        if len(data) == 0:
            # warnings.warn('SQL query get an empty')
            return data

        # fix the python <--> sql indexes
        # if atnames == 'rowID':
        if 'rowID' in columns:
            index = columns.split(',').index('rowID')
            for i in range(len(data)):
                data[i][index] -= 1

        # postporcess the output of the SQl query
        # flatten it if each els is of size 1
        if len(data[0]) == 1:
            data = [d[0] for d in data]

        return data

    def update(self, columns, values, tablename='ATOM', **kwargs):
        """Update the database with given values.

        Args:
            columns (str): names of column to update, e.g. "x,y,z".
                Avaiable columns:
                    serial, name, altLoc, resName, chainID, resSeq,
                    iCode, x, y, z, occ, temp, element, model
            values (np.ndarray): an array of values that corresponds
                        to the number of columns and atoms selected.
            kwargs: selection arguments,
                eg: name = ['CA', 'O'], chainID = ['A'],
                or no_name = ['CA', 'C'], no_chainID = ['A'].
        Examples:
            >>> values = np.array([[1.,2.,3.], [4.,5.,6.]])
            >>> db.update("x,y,z", values=values, resName='MET', name=['CA', 'CB'])
        """
        # check arguments format
        valid_colnames = self.get_colnames()

        if not isinstance(columns, str):
            raise TypeError("argument columns must be str")

        if columns != '*':
            for i in columns.split(','):
                if i not in valid_colnames:
                    raise ValueError(
                        f'Invalid column name {i}. Possible names are\n'
                        f'{self.get_colnames()}')

        # the asked keys
        keys = kwargs.keys()

        # handle the multi model cases
        if 'model' not in keys and self._nModel > 0:
            for iModel in range(self._nModel):
                kwargs['model'] = iModel
                self.update(columns, values, **kwargs)
            return

        # parse the attribute
        if ',' in columns:
            columns = columns.split(',')

        if not isinstance(columns, list):
            columns = [columns]

        # check the size
        natt = len(columns)
        nrow = len(values)
        ncol = len(values[0])

        if natt != ncol:
            raise ValueError(
                'Number of cloumns does not match between argument columns and values')

        # get the row ID of the selection
        rowID = self.get('rowID', **kwargs)
        nselect = len(rowID)

        if nselect != nrow:
            raise ValueError(
                'Number of data values incompatible with the given conditions')

        # prepare the query
        query = 'UPDATE {tablename} SET '.format(tablename=tablename)
        query = query + ', '.join(map(lambda x: x + '=?', columns))
        query = query + ' WHERE rowID=?'

        # prepare the data
        data = []
        for i, val in enumerate(values):

            tmp_data = [v for v in val]

            # here the conversion of the indexes is a bit annoying
            tmp_data += [rowID[i] + 1]

            data.append(tmp_data)

        self.c.executemany(query, data)

    def update_column(self, colname, values, index=None, tablename='ATOM'):
        """Update a single column.

        Args:
            colname (str): name of the column to update
                Avaiable columns:
                    serial, name, altLoc, resName, chainID, resSeq,
                    iCode, x, y, z, occ, temp, element, model
            values (list): new values of the column
            index (None, optional): index of the column to update (default all)

        Example:
            >>> db.update_column('x',np.random.rand(10),index=list(range(10)))
        """
        if index is None:
            data = [[v, i + 1] for i, v in enumerate(values)]
        else:
            data = [[v, ind + 1] for v, ind in zip(values, index)]

        query = 'UPDATE {tablename} SET {cn}=? WHERE rowID=?'.format(tablename=tablename,
                                                                     cn=colname)
        self.c.executemany(query, data)

    def add_column(self, colname, coltype='FLOAT', value=0, tablename='ATOM'):
        """Add an new column to the ATOM table with same value for each row.

        Args:
            colname (str): the new column name
            coltype (str): data type of the new column. Default to float.
            value (int/float/str): value for the new column.
                Default to 0.

        Examples:
            >>> db.add_column('id', coltype='str', value='positive')
        """

        query = "ALTER TABLE %s ADD COLUMN '%s' %s DEFAULT %s" % (tablename,
                                                                  colname, coltype, str(value))
        self.c.execute(query)

    def _commit(self):
        """Commit to the database."""
        self.conn.commit()
