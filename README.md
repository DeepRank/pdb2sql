# PDB2SQL

[![PyPI](https://img.shields.io/pypi/v/pdb2sql)](https://pypi.org/project/pdb2sql/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3232887.svg)](https://doi.org/10.5281/zenodo.3232887)
[![RSD](https://img.shields.io/badge/RSD-pdb2sql-red)](https://research-software.nl/software/pdb2sql)
![Build_Test](https://github.com/DeepRank/pdb2sql/workflows/Build_Test/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/pdb2sql/badge.svg)](https://coveralls.io/github/DeepRank/pdb2sql)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/36ad228df234488ab70ade6b2a80d54b)](https://www.codacy.com/gh/DeepRank/pdb2sql/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DeepRank/pdb2sql&amp;utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/pdb2sql/badge/?version=latest)](https://pdb2sql.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02077/status.svg)](https://doi.org/10.21105/joss.02077)

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
from pdb2sql import pdb2sql
pdb = pdb2sql('1AK4.pdb')
atoms = pdb.get('x,y,z',
                name = ['C','N', 'O'],
                no_resName = ['VAL','LEU'],
                chainID = 'A')
```
