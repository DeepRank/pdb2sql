# PDB2SQL

pdb2sql allows to use SQL queries to handle PDB files.
The project grew out of the developement of DeepRank and is still very much in development

At the moment two strategies are developped one using SQLite3 and the other SQLalchemy
SQLalchemy allows to have a object oriented approach but seems a bit slower.


## Installation

Clone the repository : `git clone https://github.com/DeepRank/pdb2sql`
Go in the repo and type : `pip install -e ./`
Test by going in the test folder and type : `pytest`

## pdb2sql

The following script loads the pdb file '1AK4.pdb' (must be in the same folder than the script) in a SQLite3 data base in about 0.02 seconds. You can query the data base using the pdb2sql.get(attribute,kwargs) methods

```python
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



### SQL Queries

SQL queries are quite versatile and can be used to return any attribute of the atoms with rather complex selections. As an example:

```python
xyz = db.get('x,y,z',chainID='A',resName=['VAL','LEU'],name=['CA','C','O','N'])
```

returns the positon of the CA, C, N and O atoms of all the residues 'VAL' and 'LEU' of chain A. Any other attribute can be returned (chainID, resName, name .... ) buy using it in the first argument. For example

```python
data = db.get('name,resSeq,resName',chainID='A')
```

returns the name, residue number and residue name of all the atoms in chain A.

Negative conditions can also be used to exclude some specific atoms from the selection. For example:

```python
data = db.get('name,resSeq,resName',chainID='A',no_name=['H','N'])
```

returns the name, residue number and residue name of all the atoms in chain A *except the Hydrogen and Nitrogens*. All the condition starting with "no_" are considered as negation. Therefore:

```python
data = db.get('name,resSeq,resName',chainID='A',no_resName=['VAL','LEU'])
```

will exclude the LEU and VAL residues from the selection.

### Modify the database

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
db = interface('1AK4.pdb')
contact_atoms = db.get_contact_atoms()
contact_residues = db.get_contact_residues()
```

The methods get_contact_atoms() returns here the rowID of the contact atoms. A few options are available to define the interface.

## Transform

Contains simple tranformation of the pdb such as translation and rotation. 

```python
t0 = time()
db = pdb2sql('5hvd.pdb')
print('SQL %f' %(time()-t0))

tr = np.array([1,2,3])
translation(db,tr,chainID='A')
```

