# PDB2SQL

pdb2sql allows to use SQL queries to handle PDB files. 
The project grew out of the developement of DeepRank

At the moment two strategies are developped one using SQLite3 and the other SQLalchemy
SQLalchemy allows to have a object oriented approach but seems a bit slower.

## pdb2sql

To have a feeling of how the routine works go to /pdb2sql/ and enter

```
python -i pdb2sql.py
```

This loads the pdb file '5hvd.pdb' in a SQLite3 data base in about 0.02 seconds
You can query the data base using the pdb2sql.get(attribute,kwargs) methods

For example:

```python
xyz = db.get('x,y,z',chainID='A',resName=['VAL','LEU'],name=['CA','C','O','N'])
```

returns the positon of the CA, C, N and O atoms of all the residues 'VAL' and 'LEU' of chain A.

The values of the data base can also be update with the pdb2sql.update(attribute,values,kwargs) method. For example

```python
xyz = db.get('x,y,z',chainID='A',resSeq=1)

xyz = np.array(xyz)
xyz -= np.mean(xyz)

db.update('x,y,z',xyz,chainID='A',resSeq=1)

```
Translate the residue of resSeq 1 of chain A to the center of the coordinate. 

## pdb2sqlAlchemy

To have a feeling of the alchemy interface enter 

```
python -i pdb2sqlAlchemy.py
```

Here as well you can get values from the database and update values to the data base with the methods .get() and .update(). There are some subtle difference

First as you can see it takes longer to create the data base since you need about 0.3 sec to create it. But thanks to sqlalchemy we can use negative conditions such as 

```python
# extract the xyz position of all VAL and LEU resiues of chain A but not the H atoms
xyz = db.get('x,y,z',chainID='A',resName=['VAL','LEU'],no_name=['H'])
```

Al the condition starting with "no_" are considered as negation. Hence here all the atoms are returned but the hyrdogens. We can also update the values of the data base 


```python
# put the data back 
db.update('x,y,z',xyz,chainID='A',resName=['VAL','LEU'],no_name=['H'])
```

The syntax is very similar as you can see. The main difference is the possibility to to return ATOM objects. This is achieved when no attributes are specified in the get() call

```python
atoms = db.get(chainID='A',resName=['VAL','LEU'],no_name=['H'])
```

This returns a list of ATOM object. The ATOM class is also defined in pdb2sqlAlchemy.py. We can there extract information about these atoms by calling their attributes

```python 
for at in atoms:
	print(at.name,at.x,at.y,at.z)
```

## Interface 

interface.py contains a class that subclass either pdb2sql or pdb2sql_alchemy. It allows to analyze the properties of the interface between two chains contained in the pdb file. We can extract the contact atoms, contact residues .... 


```
python -i interface.py
```

You can decide which flavor of pdb2sql to use by uncommenting/commenting the import at the top of the file

```python
#from pdb2sql import pdb2sql
from pdb2sqlAlchemy import pdb2sql_alchemy as pdb2sql
```

The definition of the class instance is pretty simple and will use either sql or alchemy

```python
db = interface('1AK4.pdb')
db.get_contact_atoms()
```

The methods get_contact_atoms() returns here the rowID of the contact atoms. A few options are available to define the interface.