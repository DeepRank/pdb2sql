======================
PDB2SQL Alchemy (beta)
======================

This module is based on `SQLAlchemy`_ that combines sql queries and object oriented programming.
Therfore :class:`~pdb2sql.pdb2sqlAlchemy.pdb2sql_alchemy` works in the same way as
:class:`~pdb2sql.pdb2sqlcore.pdb2sql` but returns arrays of objects instead of nested lists. It is however a bit slower.

.. _SQLAlchemy: https://www.sqlalchemy.org/

.. note::
    This module is not tested, use with care.

Examples
^^^^^^^^

You can get values from the database and update values to
the database with the methods .get() and .update().
The syntax is identical to the the one of pdbsqlcore:

.. ipython:: python
    :suppress:

    # change working directory to docs/
    import os
    os.chdir('..')

.. ipython:: python

    from pdb2sql.pdb2sqlAlchemy import pdb2sql_alchemy

    # create the database
    db = pdb2sql_alchemy("./pdb/dummy.pdb")

    # extract the xyz coordinates of all MET and LEU resiues of chain A but not the H atoms
    xyz = db.get('x,y,z', chainID='A', resName=['MET','LEU'], no_name=['H'])
    xyz

    # put the data back
    db.update('x,y,z', xyz, chainID='A', resName=['MET','LEU'], no_name=['H'])

The main difference is the possibility to return ATOM objects.
This is achieved when no attributes are specified in the .get() call.

.. ipython:: python

    # extract atoms
    atoms = db.get(chainID='A', resName=['MET','LEU'], no_name=['H'])

    for at in atoms:
        print(at.name,at.x,at.y,at.z)

It returns a list of :class:`~pdb2sql.pdb2sqlAlchemy.ATOM` objects.
We can there extract information about these atoms by calling their attributes.

.. automodule:: pdb2sql.pdb2sqlAlchemy
   :members:
   :undoc-members:
   :show-inheritance:
