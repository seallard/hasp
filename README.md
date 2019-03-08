# hasp
A Haskell program that builds phylogenetic trees from DNA sequences using progressive alignment and UPGMA clustering.
Developed with Ivar Mira and Alfred Nilsson as a project in a course (Program Design and Data Structures).

## Todo
* Fix buggy traceback
* Fix labeling of sequences in final tree
* Score gap opening and extension differently
* Incorporate fasta file parser
* Implement nearest neighbour clustering

## Installing and using project
To run the program, load main.hs into GHCI and run:
```
treeBuilder sequences
```
where sequences is a list of DNA strings to construct the phylogenetic tree from.
