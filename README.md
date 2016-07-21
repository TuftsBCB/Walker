# Walker

Requires modules `networkx`, `numpy`, `scikit-learn`, and `argparse`.

For a description of the Random Walk with Restart (RWR) algorithm, see the paper by Kohler et al. at [http://www.sciencedirect.com/science/article/pii/S0002929708001729](http://www.sciencedirect.com/science/article/pii/S0002929708001729).

This module was initially created to run node removal experiments with two separate graphs, but the code in matrix\_main.py can easily be used to run the standard RWR algorithm on a single graph of any sort.

## Running a random walk

The matrix\_main.py script can be used to run a random walk. The syntax looks like:

`python matrix_main.py <input_graph> <seed> [-l <low_list>] [-r <remove_nodes>]`

where the input graph is in edge list format, the seed is a list of nodes to
start the random walk at, the optional low list is a list of nodes to down-weight
for node removal experiments (as in the tissue-specific networks paper), and the
optional node removal list is a list of nodes to remove completely from the graph.

More thorough documentation to come.

## Using the module

If you use the Walker module, please cite:

Zhang H, Schaefer M, Crawford J, Kiel C, Serrano L, and Cowen LJ. "Studying Gene Prioritization in Tissue-Specific Networks: A Case Study with the BRAF Oncogene." In review. (2016)
