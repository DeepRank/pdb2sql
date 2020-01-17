---
title: 'The pdb2sql Python Package: Parsing, Manipulation and Analysis of PDB Files Using SQL Queries'
tags:
  - Python
  - Bioinformatics
  - PDB files
authors:
  - name: Nicolas Renaud
    orcid: 0000-0001-9589-2694
    affiliation: 1
  - name: Cunliang Geng
    orcid: 0000-0002-1409-8358
    affiliation: 1
affiliations:
 - name: Netherlands eScience Center, Science Park 140 1098 XG Amsterdam, the Netherlands
   index: 1
date: 13 January 2020
bibliography: paper.bib

---

# Summary

The analysis of biomolecular structures is a crucial task for a wide range of applications ranging from drug design to protein engineering. The Protein Data Bank (PDB) file format [@pdb] is the most popular format to describe biomolecular structures such as proteins and nucleic acids. In this text-based format each line represents a given atom and entails its main properties such as atom name and identifier, residue name and identifier, chain identifier, coordinates, etc. Several solutions have been developed to parse PDB files in dedicated objects that facilitate the analysis and manipulation of biomolecular structures. This is for example the case of the BioPython [@biopython,@biopdb] parser that loads PDB files in a nested dictionary whose structure mimics the hierarchical nature of the biomolecular structure. Selecting a given sub-part of the biomolecule can then be done by going through the dictionary and selecting the required atoms. Other packages, such as ``ProDy`` [@prody], ``BioJava`` [@biojava], ``MMTK`` [@mmtk] and ``MDAnalysis`` [@mdanalysis] to cite a few, also offer solution to parse PDB files. However these parsers are embedded in large applications that are sometimes difficult to integrate in new applications and are often geared toward the analysis of molecular dynamics simulations. Light-weight applications such as ``pdb-tools`` [@pdbtools] lack the capabilities to manipulate coordinates.


We present here the Python package ``pdb2sql``, which loads individual PDB files in a relational database. Among different solutions the Structured Query Language (SQL) is the most popular solution to query a given database. However SQL is its own language and domain scientists such as bioinformaticians are usually not familiar with it, which represents an important barrier for the adoption of SQL technology in bioinformatics. ``pdb2sql`` exposes complex SQL queries through simple Python methods that are intuitive for end users. Native ``SQLite3`` database can be used or ``sqlalchemy`` [@sqlalchemy] can also be leverage to obtain an object oriented approach. As such our package leverages the power of SQL queries and remove the barrier that SQL complexity represents. ``pdb2sql`` is able to parse and manipulate all atom properties of the PDB file. Besides, several advanced modules have also been built, for example to rotate or translate biomolecular structures, to characterize interface contacts, and to measure structure similarity between two protein complexes. Additional modules can easily be developed following the same scheme. As a consequence, ``pdb2sql`` is a light-weight and versatile PDB tool that is easy to extend and to integrate with new applications.

# Capabilities of ``pdb2sql``

``pdb2sql`` allows to query, manipulate and process PDB files through a series of dedicated classes. We give an overview of these features and illustrate them with snippets of code. More examples can be found in the documentation (https://pdb2sql.readthedocs.io).

## Extracting data from PDB files

``pdb2sql`` allows to simply query the database using the ``get(attr, **kwargs)`` method. The attribute ``attr`` is here a list of or a single column name of the ``SQL`` database, see Table 1 for available attributes. The keyword argument ``kwargs`` can then be used to specify a sub-selection of atoms.

Table 1. Atom attributes and assoicated definitions in ``pdb2sql``

| Attribute | Definition                                |
|-----------|-------------------------------------------|
| serial    | Atom serial number                        |
| name      | Atom name                                 |
| altLoc    | Alternate location indicator              |
| resName   | Residue name                              |
| chainID   | Chain identifier                          |
| resSeq    | Residue sequence number                   |
| iCode     | Code for insertion of residues            |
| x         | Orthogonal coordinates for X in Angstroms |
| y         | Orthogonal coordinates for Y in Angstroms |
| z         | Orthogonal coordinates for Z in Angstroms |
| occ       | Occupancy                                 |
| temp      | Temperature factor                        |
| element   | Element symbol                            |
| model     | Model serial number                       |


Every attribute name can be used to select or exclude specific atoms and multiple conditions can be easily combined. For example,

```python
from pdb2sql import pdb2sql
pdb = pdb2sql('1AK4.pdb')
atoms = pdb.get('x,y,z',
                name=['C','H'],
                resName=['VAL','LEU'],
                chainID='A')
```

this snippet extracts the coordinates of the carbon and hydrogen atoms that belong to VAL or LEU residue of the chain A.

Atoms can be excluded from the selection by appending the prefix ``no_`` to the attribute name. For example,

```python
from pdb2sql import pdb2sql
pdb = pdb2sql('1AK4.pdb')
atoms = pdb.get('name, resName',
                no_resName=['GLY', 'PHE'])
```
it extracts the atom and residue names of all atoms except those belonging to GLY or PHE residue.

## Manipulating PDB files

The data contained in the SQL database can also be modified using the ``update(attr, vals, **kwargs)`` method. The attributes and keyword arguments are identical to those in the ``get`` method. And the ``vals`` argument should contain an array of data matching the selection criteria.  For example,

```python
import numpy as np
from pdb2sql import pdb2sql
pdb = pdb2sql('1AK4.pdb')
xyz = pdb.get('x,y,z', chainID='A', resSeq=1)
xyz = np.array(xyz)
xyz -= np.mean(xyz)
pdb.update('x,y,z', xyz, chainID='A', resSeq=1)
```

this snippet first extracts the coordinates of atoms in the first residue of chain A, then translates this fragment to the origin and updates the coordinate values in the database.

``pdb2sql`` also provides a convinent class ``transform`` to easily tanslate or rotate strcutures. For example, to translate strcuture 5Å along the Y-axis,

```python
import numpy as np
from pdb2sql import pdb2sql
from pdb2sql import transform

pdb = pdb2sql('1AK4.pdb')
trans_vec = np.array([0,5,0])
transform.translation(pdb, trans_vec)
```

or to rotate strcutures 180 degrees along the X-axis,
```python
angle = np.pi
axis = (1., 0., 0.)
transform.rot_axis(pdb, axis, angle)
```

## Identifying interface

The ``interface`` class is derived from the ``pdb2sql`` class and offers functionalities to identify contact atoms or residues between two different chains with a given contact distance. It is useful for extracting and analysing the interface of e.g. protein-protein complexes.

The following example snippet returns all the atoms and all the residues of the interface of '1AK4.pdb' defined by a contact distance of 6Å.

```python
from pdb2sql import interface
pdb = interface('1AK4.pdb')
atoms = pdb.get_contact_atoms(cutoff=6.0)
res = pdb.get_contact_residues(cutoff=6.0)
```

## Computing Structure Similarity

The ``StructureSimilarity`` class allows to compute similarity measures between two protein-protein complexes. Several popular measures used to classifiy qualities of protein complex structures in the CAPRI (Critical Assessment of PRedicted Interactions) challenges [@capri] have been implemented: interface rmsd, ligand rmsd,  fraction of native contacts and DockQ[@dockq]. The approach used to compute the interface rmsd and ligand rmsd is identical to the well-known package ``ProFit`` [@profit]. All the methods required to superimpose structures have been implemented in the ``transform`` class and therefore relies on no external dependencies. The following snippet shows how to compute these measures,

```python
from pdb2sql import StructureSimilarity
sim = StructureSimilarity(decoy = '1AK4_model.pdb',
                          ref = '1AK4_xray.pdb')

irmsd = sim.compute_irmsd_fast()
lrmsd = sim.compute_lrmsd_fast()
fnat = sim.compute_fnat_fast()
dockQ = sim.compute_DockQScore(fnat, lrmsd, irmsd)
```


# Application
``psb2sql`` has been used at the Netherlands eScience center for bioinformatics projects. This is the case of ``iScore`` [@iScore] that uses graph kernels and support vector machines to rank protein-protein interface. We illustrate here the use of the pacakge by computing the interface rmsd and ligand rmsd of a series of structural models using the experimental structure as a reference. This is a common task for protein-protein docking where a large number of docked conformations are generated and have then to be compared to a ground truth to identify the best generated poses. This calculation is usuall done using the ProFit software and we therefore compare our results with those obtained with ProFit. The code do compute the similarity measure for different decoys is simple:

```python
from pdb2sql import StructureSimilarity

ref = '1AK4.pdb'
decoys = os.listdir('./decoys')
irmsd = {}

for d in decoys:
    sim = StructureSimilarity(d, ref)
    irmsd[d] = sim.compute_irmsd_fast(method='svd', izone='1AK4.izone')
```

Note that the method will compute the i-zone, i.e. the zone of the proteins that form the interface in a similar way with ProFit. This is done for the first calculations and the i-zone is then reused for the subsequent calculations. The comparison of our interface rmsd values to those given by ProFit are shown in Fig 1.

![Example figure.](sim.png)
Figure 1. Superimposed model (green) and reference (cyan) structures and comparison of interface rmsd values between pdb2sql and ProFit.

# Acknowledgements
We acknowlege contributions from Li Xue, Sonja Georgievska and Lars Ridder.


# References