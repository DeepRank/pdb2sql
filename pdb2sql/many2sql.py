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
            pdbfiles = [pdbfiles]
        self.npdb = len(pdbfiles)

        self.tablenames = tablenames
        if self.tablenames is None:
            self.tablenames = ['ATOM']
            for i in range(1, self.npdb):
                self.tablenames.append('ATOM'+str(i))

        super().__init__(pdbfiles[0], tablename=tablenames[0])

        for i in range(1, self.npdb):
            self._create_table(pdfiles[i], tablename=tablenames[i])
