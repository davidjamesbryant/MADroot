# MADroot
Efficient implementation of the MAD method for inferring the location of the root in a large phylogeny.

The MAD  method for rooting a phylogeny was first described in 

Tria, F. D. K., Landan, G., and Dagan, T. (2017). Phylogenetic rooting using minimal ancestor deviation. _Nature Ecology and Evolution_, 1, 0193.
   
This is an implementation of the O(n^2) time and linear space algorithm outlined in 

Bryant, David, and Michael Charleston. "MAD roots for large trees." _arXiv preprint arXiv:1811.03174_ (2018).

---
## Contents

*   'madRoot.cpp'  Code for the MAD root procedure, and for the simulations performed in the paper.
*   'madRoot.R'  An R script giving an (inline) version of the madRoot algorithm
*   'Makefile'     Make file.
*   'smallExample.txt'  Small (four taxa) example
*    'largeExample.txt' Large (3082 taxa) example
Note that this code makes use of the PhyLib library (here included as a submodule).

---
## Downloading

In a terminal window, execute
    git clone --recursive git@github.com:davidjamesbryant/MADroot.git
    
---
## Building

### madRoot executable
Change to the directory 'MADroot'. Make sure that 'madRoot.cpp', 'Makefile' and the 'Phylib' directory are in the same folder, then type 'make'.

To run the algorithm on the tree in the tree file <filename> call
    madRoot <treefile>

To repeat the simulation/benchmark performed in Bryant and Charleston, run
    madRoot -SIM

### madRoot (inline) R.

in R, make MADroot your working directory. Install the packages 'Rcpp' and 'inline'. When you execute 'madRoot'R' it will create a new function in R called 'madRoot'. You can run 'madRoot' using
    madRoot(_treeString_)
where _treeString_ is a string (character vector) containing one or more trees in NEWICK notation, separated by semicolons.  
