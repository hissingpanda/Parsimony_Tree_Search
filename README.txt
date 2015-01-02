Created by: 
Paul Lanctot
lanct020

TOPIC III: Parsimony Tree Search Algorithm
-------------------------------------------

Matlab files in this directory:
sankoffAmino.m
treeParseAmino.m
parsimonyTreeSearch.m
nni.m
parsimonyTreeSearch2.m
nni2.m

Extra files in this directory:
Final_Report_lanctot.pdf
Extras/
	Phylogenetic_trees(.)_seed.jpg
	Phylogenetic_trees(.)_full.jpg
	Phylogenetic_trees(.)_seed.tree
	Phylogenetic_trees(.)_full.tree
	scores_seed.txt
	scores_full.txt	

Matlab files
-------------

-----------------------
sankoffAmino.m
----------------------
Function sankoffAmino that performs Sankoff's algorithm on two phylogenetic nodes
         Input: Left and right node of a branch
        Output: The parsimony values for each nucleotide at the node 

You can call the function in Matlab by:
sankoffAmino(left_node_of_tree_branch, right_node_of_tree_branch)


--------------------------
treeParseAmino.m
-------------------------
Function treeParseAmino that performs Sankoff's algorithm on a phylogenetic tree.
     Input: phylogenetic tree, multiple sequence alignment
    Output: the parsimony score of the input tree, the sequence alignment of each internal node/branch

You can call the function in Matlab by:
treeParseAmino(phylogenetic_tree, alignment)


--------------------------
nni.m
-------------------------
Function nni that searches for the best possible phylogenetic tree, using
    the nearest-neighbor interchange algorithm
         Input: phylogenetic tree, parsimony score, multiple sequence alignment
        Output: phylogenetic tree, parsimony score, value to see if score changed

You can call the function in Matlab by:
nni(phylogenetic_tree, parsimony, alignment)


--------------------------
parsimonyTreeSearch.m (Main file: run this to execute all the tasks)
-------------------------
Function parsimonyTreeSearch that performs three tasks:

        1) Builds initial phylogenic tree using pairwise distances of the
            input multiple sequence alignment
        2) Computes the tree's parsimony score
        3) Searches for a better tree (and generates it) using the parsimony score

         Input: multiple sequence alignment
        Output: phylogenetic tree, parsimony score for the input tree

You can call the function in Matlab by:
parsimonyTreeSearch('PF02171/PF02171_seed.fasta')

Note: In order to run this on the full fasta file, need to use the commented out seqpdist(seqs) line instead of the current version.
	Otherwise a matrix indexing error will occur due to the pairwise-deletions.

---------------------
parsimonyTreeSearch2.m
nni2.m
------------------

These are cut-down versions of the aformentioned algorithms, used to produce the full PIWI family in a faster, yet still accurate, way.

They can be run in Matlab similarly to their alternative functions.

---------
Extras/
--------
scores_seed.txt: Parsimony scores at different iterations of the generated tree (using the seeds file).
scores_full.txt: Parsimony scores at different iterations of the generated tree (using the full file).
Phylogenetic_trees(.)_seed.jpg: Several image files that show different iterations of the generated tree as the program was ran (shown in report).
Phylogenetic_trees(.)_full.jpg: Several image files that show different iterations of the generated tree as the program was ran (shown in report).
Phylogenetic_trees(.)_seed.tree: Tree files of the images that can be used in matlab.
Phylogenetic_trees(.)_full.tree: Tree files of the images that can be used in matlab.
Final_Report_Lanctot.pdf: 10-page project report with brief introduction of the implementation and details of your results and answers to the questions.
