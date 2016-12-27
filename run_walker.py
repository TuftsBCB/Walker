"""
Main script for running tissue-specific graph walk experiments, to convergence.

"""
import sys
import argparse
from walker import Walker

def generate_seed_list(seed_file):
    """ Read seed file into a list. """
    seed_list = []

    try:
        fp = open(seed_file, "r")
    except IOError:
        sys.exit("Error opening file {}".format(seed_file))

    for line in fp.readlines():
        info = line.rstrip().split()
        if len(info) > 1:
            seed_list.append(info[1])
        else:
            seed_list.append(info[0])

    fp.close()
    return seed_list

def get_node_list(node_file):
    node_list = []
    try:
        fp = open(node_file, 'r')
    except IOError:
        sys.exit('Could not open file: {}'.format(node_file))

    # read the first (i.e. largest) connected component
    cur_line = fp.readline()
    while cur_line and not cur_line.isspace():
        if cur_line:
            node_list.append(cur_line.rstrip())
        cur_line = fp.readline()

    fp.close()
    return node_list

def main(argv):

    # set up argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('input_graph', help='Original graph input file, in\
                                             edge list format')
    parser.add_argument('seed', help='Seed file, to pull start nodes from')
    parser.add_argument('-e', '--restart_prob', type=float, default=0.7,
                        help='Restart probability for random walk')
    parser.add_argument('-l', '--low_list', nargs='?', default=None,
                        help='<Optional> List of genes expressed and\
                              unexpressed in the current tissue, if applicable')
    parser.add_argument('-n', '--node_list', nargs='?', default=None,
                        help='<Optional> Order of output probs')
    parser.add_argument('-o', '--original_graph_prob', type=float, default=0.1,
                        help='Probability of walking on the original (non-\
                              tissue specific) graph, if applicable')
    parser.add_argument('-r', '--remove', nargs='+',
                        help='<Optional> Nodes to remove from the graph, if any')
    opts = parser.parse_args()

    seed_list = generate_seed_list(opts.seed)
    node_list = get_node_list(opts.node_list) if opts.node_list else []

    # filter nodes we want to remove out of the starting seed, if any
    remove_list = opts.remove if opts.remove else []
    if remove_list:
        seed_list = [s for s in seed_list if s not in remove_list]

    # run the experiments, and write a rank list to stdout
    wk = Walker(opts.input_graph, opts.low_list, remove_list)
    wk.run_exp(seed_list, opts.restart_prob,
               opts.original_graph_prob, node_list)


if __name__ == '__main__':
    main(sys.argv)

