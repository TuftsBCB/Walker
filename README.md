# Walker

[![DOI](https://zenodo.org/badge/63801061.svg)](https://zenodo.org/badge/latestdoi/63801061)

Requires modules `networkx`, `numpy`, `scikit-learn`, and `argparse`.

For a description of the Random Walk with Restart (RWR) algorithm, which
this module implements, see the paper by Kohler et al. at
[http://www.sciencedirect.com/science/article/pii/S0002929708001729](http://www.sciencedirect.com/science/article/pii/S0002929708001729).

## Overview

This module can be used to run two types of experiments:

- A standard random walk with restart from a set of seed nodes, as in the
  Kohler et al. paper referenced above.
- A random walk with restart, from a set of seed nodes, on a "tissue-specific"
  network. The network is defined by a "low list" of nodes (i.e. genes) that
  are not expressed in the tissue of interest. This is described in more
  detail in our paper, which is currently in review.

Examples of both experiments are described in more detail below.

## Running a random walk

The run\_walker.py script can be used to run a random walk. The syntax looks like:

`python run_walker.py <input_graph> <seed> [-l <low_list>] [-r <remove_nodes>]`

where the input graph is in edge list format, the seed is a list of nodes to
start the random walk at, the optional low list is a list of nodes to down-weight
for node removal experiments, and the optional node removal list is a list of nodes
to remove completely from the network.

The script will write a tab-separated list of nodes and probabilities to stdout,
where the probability number represents the probability that a random walk
starting at the seed nodes will terminate at the given node.

For more detail about the expected arguments, run `python run_walker.py -h`.

## Examples

To help you get up and running, a few simple examples are included in the `testdata`
folder. To run a standard random walk experiment on a simple example network, run
this command:

`python run_walker.py testdata/test_network.ppi testdata/test_seed.txt`

Or, to run a "tissue-specific" random walk experiment using the same
simple example network, try:

`python run_walker.py testdata/test_network.ppi testdata/test_seed.txt -l testdata/test_low_list.txt`

## Using the module

If you use the Walker module, please cite:

Zhang H, Schaefer M, Crawford J, Kiel C, Serrano L, and Cowen LJ. "Studying Gene Prioritization in Tissue-Specific Networks: A Case Study with the BRAF Oncogene." In review. (2016)
