"""
PDB2SQL

A package to leverage SQL queries to parse, manipulate and process PDB files.

Provides:
1. a powerful pdb2sql object to convert PDB data in SQL database
2. strcuture transformation functions (rotations, translations...)
3. useful capablities to
    - calculate structure interface (contact atoms and residues)
    - calculate structure similarity (iRMSD, lRMSD, FNAT, DockQ...)

Reference:
    Online tutorial and documentation: https://pdb2sql.readthedocs.io

Example:
    `pdb2sql` easily allows to load a PDB file in an object. Once loaded,
    the data can be parsed using SQL queries. To facilitate the adoption of
    the tool simple methods have been developped to wrap the SQL queries in
    simple methods.

    For example obtaining the positions of all carbon, nitrogen and oxygen
    atoms of chain A from all residues but VAL and LEU, one can use :

    >>> from pdb2sql import pdb2sql
    >>> pdb = pdb2sql('1AK4.pdb')
    >>> atoms = pdb.get('x,y,z',
    ...                 name = ['C','N', 'O'],
    ...                 no_resName = ['VAL','LEU'],
    ...                 chainID = 'A')


Available modules:
    pdb2sql
        Core `pdb2sql` object
    many2sql
        Core `many2sql` object
    interface
        Core `interface` object
    StructureSimilarity
        Tools to compute structure similarities between two structures.
    transform
        Tools to do structure transformation
    align
        Tools to do structure alignment
    superpose
        Tools to do structure superposition

Utilities:
    fetch
        download PDB file from PDB website https://www.rcsb.org/.
"""

from .pdb2sqlcore import pdb2sql
from .many2sql import many2sql
from .interface import interface
from .StructureSimilarity import StructureSimilarity
from . import transform
from .utils import fetch
from .align import align, align_interface
from .superpose import superpose

from .__version__ import __version__

# remove unnecesary modules
del pdb2sql_base
del pdb2sqlcore
del utils
