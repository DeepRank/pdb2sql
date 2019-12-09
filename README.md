# PDB2SQL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3232888.svg)](https://doi.org/10.5281/zenodo.3232888)
[![Travis](https://secure.travis-ci.org/DeepRank/pdb2sql.svg?branch=master)](https://travis-ci.org/DeepRank/pdb2sql)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/pdb2sql/badge.svg)](https://coveralls.io/github/DeepRank/pdb2sql)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7cce335b21cc49d6bb3ff9878d6ac904)](https://www.codacy.com/manual/CunliangGeng/pdb2sql?utm_source=github.com&utm_medium=referral&utm_content=DeepRank/pdb2sql&utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/pdb2sql/badge/?version=latest)](https://pdb2sql.readthedocs.io/en/latest/?badge=latest)

PDB2SQL is a Python package that allows to use SQL queries to handle PDB files.
This project grew out of the developement of DeepRank.

-   Source code: <https://github.com/DeepRank/pdb2sql>
-   Documentation: <https://pdb2sql.readthedocs.io>

It provides:

-   a powerful `pdb2sql` object to manipulate PDB data in SQL database
-   strcuture transformation functions (rotations, translations...)
-   useful capablities to
    -   calculate structure interface (contact atoms and residues)
    -   calculate structure similarity (iRMSD, lRMSD, FNAT, DockQ...)
