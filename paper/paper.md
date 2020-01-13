---
title: 'The pdb2sql Python Module : Parsing, Manipulation and Analysis of PDB Files Using SQL Queries'
tags:
  - Python
  - Bioinformatics
  - PDB files
authors:
  - name: Nicolas Renaud
    orcid: 0000-0001-9589-2694
    affiliation: 1 
  - name: Cunliang Geng
    orcid:
    affiliation: 1
affiliations:
 - name: Netherlands eScience Center, Science Park 140 1098 XG Amsterdam, the Netherlands
   index: 1
date: 13 January 2020
bibliography: paper.bib

---

# Summary


The analysis of protein structure is a crucial task for a wide range of applications ranging from drug design to enzyme catalyziz. The Protein Data Bank (PDB) file format is the most popular format to describe protein structure. In this text-based format each line represents a given atom and entails its main properties such as atom type, positions, residues etc .... Several solutions have been developed to parse  PDF files in dedicated objects that facilitate the analysis and manipulation of protein structure. This is for example the case of the BioPython [@biopython] parser that loads PDB files in a nested dictionary whose structure mimic the hierarchical nature of the protein structure. Selecting a given sub-part of the protein can then be done by going through the dictionaries and selecting the required atoms.

Relational database are extensively used to organize data and to query data banks with impressive performance. Among different solutions the Structured Query Language is the most popular solution to query a given database. However SQL is its own language and domain scientists such as bioinformatician are not exposed to it. This lack of exposure represents an important barrier for the adoption of SQL technology in bioinformatics.

We present here Python package called ``pdb2sql`` that laods individual PDB files in a SQL database and expose complex SQL queries through simple methods that are intuitive for end users. As such our library leverage the power of SQL queries and remove the barrier that SQL complexity represents. Native ``SQLlite3`` database can be used or ``sqlalchemy`` [@sqlalchemy] can be leverage to obtain an object oriented approach. Additional modules have been build on top of this technology for example to characterize protein-protein interfaces and to mesure similarity structure between two proteins. Additional modules can easily be developed following the same scheme.


## Extracting data from a PDB files

``pdb2sql`` allows to simply query the database using the ``.get(attr,**kwargs)`` method. The attribute is here a list of or a single column name of the ``SQL`` database. Choices are: 
  * serial : index of the atom
  * name : atom name
  * altLoc : ,
  * resName : name of the residue 
  * chainID : name of the chain
  * resSeq : residue index
  * x,y,z : coordinate of the atom
  * occ : occupation number
  * temp : temperature
  * model : model number


The keyword argument can then be used to specify a sub selection of atoms. Every attribute name can be used to select or exclude specific atoms and multiple conditions can be easily combined. For example the following snippets

```python
from pdb2sql.pdb2sqlcore import pdb2sql
pdb = pdb2sql('1AK4.pdb')
pdb.get('x,y,z',name=['C','H'], 
         resName=['VAL','LEU'], chainID='A')
```

extracts the positions of the carbon and hydrogen atoms that belong to a VAL or LEU residue of the chain A. Atoms can be excluded from the selection by appending the prefix ``no_`` to the attribute name.For example the following snippet extracts all  atom and residue names of all atoms except those belonging to a GLY or a PHE residue.

```python
from pdb2sql.pdb2sqlcore import pdb2sql
pdb = pdb2sql('1AK4.pdb')
pdb.get('name, resName', no_resName=['GLY', 'PHE'])
```

## Manipulating PDB Files

The data contained in the sql data base can also be modified using the ``.update(attr,vals,**kwargs)`` method. The attributes and keyword arguments are identical to the ``.get(attr,vals,**kwargs)`` method and the ``vals`` argument should contain an array of data matching the selection criteria.  For example the following code first extract the positions of atoms in the first residue of chain A, then translate this fragment to the origin and update the coordinate values in the database.

```python
import numpy as np
from pdb2sql.pdb2sqlcore import pdb2sql
pdb = pdb2sql('1AK4.pdb')
xyz = pdb.get('x,y,z',chainID='A',resSeq=1)
xyz = np.array(xyz)
xyz -= np.mean(xyz)
pdb.update('x,y,z',xyz,chainID='A',resSeq=1)
```

# Identifying protein-protein interface

The ``interface()`` class is derived from the ``pdb2sql()`` class and offers functionalities to identify contact atoms of the different chains composing the structure. Contact atoms are defined as atom of a given chain with an atom of a different chain located less than a given distance away. For example the following snippet returns all the atoms and all the residue of the interface of '1AK4.pdb' defined by a contact distance of 6.0 Angstrom.

```python
from pdb2sql.interface import interface
pdb = interface('1AK4.pdb')
atoms = pdb.get_contact_atoms(cutoff=6.0)
res = pdb.get_contact_residues(cutoff=6.0)
```

## Computing Structure Similarity

The ``StructureSimilarity()`` class allows to compute similarity measure between different structures. Several measures have been implemented : interface-rmsd, ligand-rmsd, Fnat and DockQ. The approach used to compute the i-rmsd and l-rmsd is identical to the well known package ``ProFit`` [@profit]. In case of the i-rmsd the longest chains of the two structure are superposed using translation and rotation. The RMSD of the atoms belonging to the remaining chains is then evaluated. All the methods required to superposed the data have been implemented in the ``transform.py`` file and therefore relies on no external dependencies.

```python
from pdb2sql.StructureSimilarity import StructureSimilarity
sim = StructureSimilarity('1AK4w_101.pdb','1AK4.pdb')

v = sim.compute_irmsd_fast(method='svd',
            izone='1AK4.izone')
        
w = sim.compute_irmsd_pdb2sql(method='svd',
            izone='1AK4.izone')
```

# Example
``psb2sql`` has been used at the Netherlands eScience center for bioinformatics projects. This is the case of ``iScore`` [@iScore] that uses graph kernels and support vector machines to rank protein-protein interface. We illustrate here the use of the library by computing the i-rmsd and l-rmsd of a series of decoy using the experimental structure as a reference. This is a common task for protein-protein docking simulations where a large number of artificial conformations are generated and have then to be compared to a ground truth to identify the best generated poses. This calculation is usuall done using the ProFit software and we therefore compare our results with thos obtained with ProFit. The code do compute the similarity measure for different decoys is simply :

```python
from pdb2sql.StructureSimilarity import StructureSimilarity

ref = '1AK4.pdb'
decoys = os.listdir('./decoys')
irmsd = {}

for d in decoys:
    sim = StructureSimilarity(d,ref)
    irmsd[d] = sim.compute_irmsd_fast(method='svd',
            izone='1AK4.izone')
```

Note that the method will compute the i-zone, i.e. the zone of the proteins that form the interface in a similar way than ProFit. This is done for the first calculations and the zone file is then reused for the subsequent calculations. The comparison of our i-rmsd value to those given by ProFit are show in Fig.  as well as the superimposed structure of a given decoy.

![Example figure.](sim.png)

# Acknowledgements
We aknowledge contribution from Li Xue, Sonja Georgievska and Lars Ridder.


# References