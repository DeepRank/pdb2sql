# Change Log


All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## 0.5.1
- Updated `compute_CapriClass` conditions

## 0.5.0
- Added atom name selection for lrmsd calculation
- Removed hardcoded chainIDs and added chainID selection in `StructureSimilarity`

## 0.4.0
- Added `many2sql` to support read multiple PDB files
- Added support for `Path` objects of input PDB files
- Added support for `help(pdb2sql)`
- Updated assignment of chain IDs in `StructureSimilarity`

## 0.3.0
- Added `align` to superpose a structure to a specific axis or plane
- Added `superpose` to superpose two structures based on selections
- Added `fetch` to download PDB file with given PDB ID
- Updated `interface` object to take `pdb2sql` object as input