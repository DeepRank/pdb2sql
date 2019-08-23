# PDB2SQL

pdb2sql allows to use SQL queries to handle PDB files.
The project grew out of the developement of DeepRank and is still very much in development.

At the moment two strategies are developped one using SQLite3 and the other SQLalchemy.
SQLalchemy allows to have a object oriented approach but seems a bit slower.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3232888.svg)](https://doi.org/10.5281/zenodo.3232888)


[![Build Status](https://secure.travis-ci.org/DeepRank/pdb2sql.svg?branch=master)](https://travis-ci.org/DeepRank/pdb2sql)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/pdb2sql/badge.svg?branch=master)](https://coveralls.io/github/pdb2sql/iScore?branch=master)



## Installation
<!-- 
  1. Clone the repository : `git clone https://github.com/DeepRank/pdb2sql`

  2. Go in the repo and type : `pip install -e ./`

  3. Test by going in the test folder and type : `pytest` -->

`pip install pdb2sql` 

## pdb2sql

The following script loads the pdb file '1AK4.pdb' (must be in the same folder than the script) in a SQLite3 data base in about 0.02 seconds. You can query the data base using the ```pdb2sql.get(attribute,**kwargs)``` method.

```python
from pdb2sql.pdb2sqlcore import pdb2sql

#create the database
db = pdb2sql('1AK4.pdb')
print('SQL %f' %(time()-t0))

# get the xyz of all the atoms
xyz = db.get('x,y,z',model=0)

# get the xyz of all the CA, C, O, N atoms of all VAL and LEU residues of chain A
xyz = db.get('x,y,z',chainID='A',resName=['VAL','LEU'],name=['CA','C','O','N'])

# move the resiude 1 of chain A
xyz = db.get('x,y,z',chainID='A',resSeq=1)
xyz = np.array(xyz)
xyz -= np.mean(xyz)
db.update('x,y,z',xyz,chainID='A',resSeq=1)

```



#### SQL Queries

SQL queries are quite versatile and can be used to return any attribute of the atoms with rather complex selections. As an example:

```python
xyz = db.get('x,y,z',chainID='A',resName=['VAL','LEU'],name=['CA','C','O','N'])
```

returns the positon of the CA, C, N and O atoms of all the residues 'VAL' and 'LEU' of chain A. Any other attribute can be returned (chainID, resName, name .... ) buy using it in the first argument. For example

```python
data = db.get('name,resSeq,resName',chainID='A')
```
returns the name, residue number and residue name of all the atoms in chain A.

#### Negative conditions

Negative conditions can also be used to exclude some specific atoms from the selection. For example:

```python
data = db.get('name,resSeq,resName',chainID='A',no_name=['H','N'])
```

returns the name, residue number and residue name of all the atoms in chain A **except the Hydrogen and Nitrogens**. All the condition starting with ```no_``` are considered as negation. Therefore:

```python
data = db.get('name,resSeq,resName',chainID='A',no_resName=['VAL','LEU'])
```

will exclude the LEU and VAL residues from the selection.

#### Modify the database

The values of the data base can also be update with the pdb2sql.update(attribute,values,kwargs) method. For example

```python
xyz = db.get('x,y,z',chainID='A',resSeq=1)
xyz = np.array(xyz)
xyz -= np.mean(xyz)
db.update('x,y,z',xyz,chainID='A',resSeq=1)
```

Translate the residue of resSeq 1 of chain A to the center of the coordinate. Note that a dedicated module called transform.py can handle translation,rottion, etc of xyz coordinates

## pdb2sqlAlchemy

SQLalchemy combine sql queries and object oriented programming. Therfore pdb2sqlAlchemy works in the same way that pdb2sqlcore but returns arrays of objects instead of nested lists. It is however a bit slower.

```python
from pdb2sql.pdb2sqlAlchemy import pdb2sql_alchemy

#create the database
db = pdb2sql_alchemy('1AK4.pdb')

# extract the xyz position of all VAL and LEU resiues of chain A but not the H atoms
xyz = db.get('x,y,z',model=0) #chainID='A',resName=['VAL','LEU'],no_name=['H'])

# put the data back
db.update('x,y,z',xyz)

# extract atoms
atoms = db.get(chainID='A',resName=['VAL','LEU'],no_name=['H'])

for at in atoms:
	print(at.name,at.x,at.y,at.z)
```

Here as well you can get values from the database and update values to the data base with the methods .get() and .update(). The syntax is identical to the the one of pdbsqlcore:


```python
# extract the xyz position of all VAL and LEU resiues of chain A but not the H atoms
xyz = db.get('x,y,z',chainID='A',resName=['VAL','LEU'],no_name=['H'])

# put the data back
db.update('x,y,z',xyz,chainID='A',resName=['VAL','LEU'],no_name=['H'])
```

#### Return ATOM objects

The main difference is the possibility to to return ATOM objects. This is achieved when no attributes are specified in the .get() call

```python
atoms = db.get(chainID='A',resName=['VAL','LEU'],no_name=['H'])
```

This returns a list of ATOM object. The ATOM class is also defined in pdb2sqlAlchemy.py. We can there extract information about these atoms by calling their attributes

```python
for at in atoms:
	print(at.name,at.x,at.y,at.z)
```

## Interface

The module interface.py contains a class that subclass pdb2sqlcore (Test for pdb2sqlAlchemy not doneyet). It allows to analyze the properties of the interface between two chains contained in the pdb file. The class allows to easily extract the contact atoms and contact residues of the conformation.

```python
from pdb2sql.interface import interface

db = interface('1AK4.pdb')
contact_atoms = db.get_contact_atoms()
contact_residues = db.get_contact_residues()
```

The methods get_contact_atoms() returns here the rowID of the contact atoms. A few options are available to define the interface.

## Structure Similarity

The StructureSimilarity module allows to computeL `irmsd, lrmsd, Fnat` and `dockQ` score of  given conformation with respect to its native. The native can be any other conformations as long as the sequences are aligned.

```python
from pdb2sql.StructureSimilarity import StructureSimilarity

# create the class instance
sim = StructureSimilarity('1AK4_300w.pdb','1AK4.pdb')

# compute the irmsd with the two different methods
irmsd_fast = sim.compute_irmsd_fast(method='svd',izone='1AK4.izone')
irmsd = sim.compute_irmsd_pdb2sql(method='svd',izone='1AK4.izone')

# compute the lrmsd with the two different methods
lrmsd_fast = sim.compute_lrmsd_fast(method='svd',lzone='1AK4.lzone',check=True)
lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None,method='svd')

# compute the Fnat with the two different methods
Fnat_fast = sim.compute_Fnat_fast(ref_pairs='1AK4.ref_pairs')
Fnat = sim.compute_Fnat_pdb2sql()

# compute the DOCKQ
dockQ = sim.compute_DockQScore(Fnat_fast,lrmsd_fast,irmsd_fast)
```

As you can see two methods are possible for the calculation of each quantity. We recommend using the **fast** that is faster and better tested.
