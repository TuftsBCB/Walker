"""
Main script for running tissue-specific graph walk experiments, to convergence.

"""
import sys
import argparse
from multi_matrix import MultiMatrix

def generate_seed_list(seed_file):
    """ Read seed file into a list. """
    seed_list = []

    try:
        fp = open(seed_file, "r")
    except IOError:
        sys.exit("Error opening file {}".format(seed_file))

    for line in fp.readlines():
        info = line.split()
        seed_list.append(info[1])

    fp.close()
    return seed_list


def main(argv):

    # set up argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("input_graph", help="Original graph input file, in\
                                             edge list format")
    parser.add_argument("seed", help="Seed file, to pull start nodes from")
    parser.add_argument("-l", "--low_list", nargs="?", default=None,
                        help="<Optional> List of genes expressed and\
                              unexpressed in the current tissue, if applicable")
    parser.add_argument("-r", "--remove", nargs="+",
                        help="<Optional> Nodes to remove from the graph, if any")
    opts = parser.parse_args()

    restart_prob = 0.7
    original_graph_prob = 0.1

    seed_list = generate_seed_list(opts.seed)

    # filter nodes we want to remove out of the starting seed, if any
    remove_list = opts.remove if opts.remove else []
    if remove_list:
        seed_list = [s for s in seed_list if s not in remove_list]

    # run the experiments, and write a rank list to stdout
    mm = MultiMatrix(opts.input_graph, opts.low_list, remove_list)
    mm.run_exp(seed_list, restart_prob, original_graph_prob)


if __name__ == '__main__':
    main(sys.argv)
