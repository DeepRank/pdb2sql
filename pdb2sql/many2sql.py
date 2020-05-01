import sqlite3
import warnings
import subprocess as sp
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

from .pdb2sqlcore import pdb2sql


class many2sql(pdb2sql):

    def __init__(self, pdbfiles, tablenames=None):
        """Create a sql dabase containing multiple pdbs."""

        if not isinstance(pdbfiles, list):
            raise ValueError('pdbfiles must be a list')

        self.npdb = len(pdbfiles)

        self.tablenames = tablenames
        if self.tablenames is None:
            self.tablenames = ['ATOM']
            for i in range(1, self.npdb):
                self.tablenames.append('ATOM'+str(i))

        super().__init__(pdbfiles[0], tablename=self.tablenames[0])

        for i in range(1, self.npdb):
            self._create_table(
                pdbfiles[i], tablename=self.tablenames[i])

    def __call__(self, **kwargs):
        """Return a class instance containing the selection of each structure

        Returns:
            many2sql: class instance containing the selection of each structure
        """

        names = self._get_table_names()

        first = True
        for n in names:
            pdb_data = self.sql2pdb(tablename=n, **kwargs)
            if first:
                new_db = many2sql([pdb_data], tablenames=[n])
                first = False
            else:
                new_db._create_table(pdb_data, tablename=n)

        return new_db

    def intersect(self, match=['name', 'resname', 'resSeq', 'chainID']):
        """Returns a many2sql instance containing the common part of all the structures.

        Args:
            match (list, optional): column name that must match in the intersection.
                                    Defaults to ['name', 'resname', 'resSeq', 'chainID'].

        Returns:
            many2sql: a class instance containing the tables of the matchin structure
        """

        all_data = self.get_intersection('*', match=match)
        all_names = self._get_table_names()

        first = True
        for name, data in zip(all_names, all_data):
            if first:
                new_db = many2sql(
                    [self.data2pdb(data)], tablenames=[name])
                first = False
            else:
                new_db._create_table(
                    self.data2pdb(data), tablename=name)

        return new_db

    def get_all(self, columns, **kwargs):
        """Returns the data from the selection of all table in the instance

        Args:
            columns (str): column tame to return

        Returns:
            list: data per structure
        """

        names = self._get_table_names()
        data = []
        for n in names:
            data.append(self.get(columns, tablename=n, **kwargs))
        return data

    def get_intersection(self, column, match=['name', 'resname', 'resSeq', 'chainID']):
        """Return the data of the interection

        Args:
            column (str): column table to return
            match (list, optional): column name that must match in the intersection.
                                    Defaults to ['name', 'resname', 'resSeq', 'chainID'].

        Returns:
            list: data per structure
        """
        names = self._get_table_names()
        ntable = len(names)
        select = "select "

        # column names
        if column == '*':
            column_list = list(self.col.keys())
        else:
            column_list = column.split(',')
        ncol = len(column_list)

        # fields to select
        fields = ''
        for n in names:
            for c in column.split(','):
                fields += n+'.'+c+', '
        fields = fields[:-2]+' '

        # join the table
        from_join = 'from ' + ' INNER JOIN '.join(names) + ' '

        # conditions
        cond = 'on '
        for attr in match:
            for i1 in range(ntable-1):
                table1 = names[i1]
                for i2 in range(i1+1, ntable):
                    table2 = names[i2]
                    cond += table1+'.'+attr+'='+table2+'.'+attr+' and '
        cond = cond[:-5]+';'
        query = select+fields+from_join+cond
        raw_data = self.conn.execute(query)

        data = []
        for i in range(ntable):
            data.append([])

        for x in raw_data:
            for it in range(ntable):
                s, e = it*ncol, (it+1)*ncol
                data[it].append(list(x[s:e]))
        return data
