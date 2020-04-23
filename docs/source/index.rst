*************************************
pdb2sql: Processing PDB data with SQL
*************************************

`pdb2sql`_ is a Python package that allows to use SQL queries to handle `PDB`_ files.
Currently, only 'ATOM' data is parsed, and other items of PDB, e.g. HETATM, are ignored.

Installation:
    ``pip install pdb2sql``

.. _pdb2sql: https://github.com/DeepRank/pdb2sql
.. _PDB: https://www.rcsb.org/

========
Tutorial
========

.. toctree::
    :maxdepth: 1

    10 minutes to pdb2sql <tutorial.rst>

==========
Python API
==========

.. toctree::
    :maxdepth: 1

    PDB2SQL <pdb2sql.pdb2sqlcore>
    Interface <pdb2sql.interface>
    Structure Similarity <pdb2sql.StructureSimilarity>
    Structure Transformation <pdb2sql.transform>
    Utilities <pdb2sql.utils>

.. :caption: Python API
..
    ==================
    Indices and tables
    ==================

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
