# PDB2SQL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3232888.svg)](https://doi.org/10.5281/zenodo.3232888)
[![Travis](https://secure.travis-ci.org/DeepRank/pdb2sql.svg?branch=master)](https://travis-ci.org/DeepRank/pdb2sql)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/pdb2sql/badge.svg)](https://coveralls.io/github/DeepRank/pdb2sql)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7cce335b21cc49d6bb3ff9878d6ac904)](https://www.codacy.com/manual/CunliangGeng/pdb2sql?utm_source=github.com&utm_medium=referral&utm_content=DeepRank/pdb2sql&utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/pdb2sql/badge/?version=latest)](https://pdb2sql.readthedocs.io/en/latest/?badge=latest)

`pdb2sql` is a Python package that leverage SQL queries to parse, manipulate and process PDB files. It provides:

-   a powerful `pdb2sql` object to convert PDB data in SQL database
-   strcuture transformation functions (rotations, translations...)
-   useful capablities to
    -   calculate structure interface (contact atoms and residues)
    -   calculate structure similarity (iRMSD, lRMSD, FNAT, DockQ...)

## Installation

```
pip install pdb2sql
```

## Documentation
The documentation of the package alongside small tutorial can be found at :
-  <https://pdb2sql.readthedocs.io>

## Quick Example

`pdb2sql` easily allows to load a PDB file in an object. Once loaded, the data can be parsed using SQL queries. To facilitate the adoption of the tool simple methods have been developped to wrap the SQL queries in simple methods. For example obtaining the positions of all carbon, nitrogen and oxygen atoms of chain A from all residues but VAL and LEU, one can use :

```python
from pdb2sql.pdb2sqlcore import pdb2sql
pdb = pdb2sql('1AK4.pdb')
atoms = pdb.get('x,y,z', 
                name = ['C','N', 'O'], 
                no_resName = ['VAL','LEU'],
                chainID = 'A')
```
