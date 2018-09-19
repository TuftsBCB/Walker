import os
import unittest
import subprocess

class ReorderTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(ReorderTest, self).__init__(*args, **kwargs)
        self.rwr_root = '/Users/jake/Documents/research/continuous_mr/giant/'
        self.rwr_script = os.path.join(self.rwr_root, 'Walker/run_walker.py')
        self.small_network = os.path.join(self.rwr_root,
                                          'Walker/testdata/test_network.ppi')
        self.large_network = os.path.join(self.rwr_root,
                                          'networks/protein.links.human.experiment.txt')


    # @unittest.skip('ignore')
    def test_reorder_small(self):
        """ Test reorder of seed nodes on small test network """
        import itertools

        num_seed_nodes = 3
        nodes = self._get_nodes(self.small_network)
        node_list = sorted(list(nodes))[:num_seed_nodes]
        results = self._run_seed_permutations(node_list,
                                              self.small_network)

        ixs = list(range(len(results)))
        for i1, i2 in itertools.combinations(ixs, 2):
            assert (sorted(results[i1], key=lambda x: x[1]) ==
                    sorted(results[i2], key=lambda x: x[1]))


    # @unittest.skip('ignore')
    def test_reorder_large(self):
        """ Test reorder of seed on large network (e.g. PPI) """
        import itertools

        num_seed_nodes = 3
        nodes = self._get_nodes(self.large_network)
        node_list = list(nodes)[:num_seed_nodes]
        results = self._run_seed_permutations(node_list,
                                              self.large_network,
                                              delim=' ')

        ixs = list(range(len(results)))
        for i1, i2 in itertools.combinations(ixs, 2):
            assert (sorted(results[i1], key=lambda x: x[1]) ==
                    sorted(results[i2], key=lambda x: x[1]))


    def _run_seed_permutations(self, node_list, network_file, delim='\t'):
        import itertools, tempfile
        results = []
        for nl in itertools.permutations(node_list):
            f = tempfile.NamedTemporaryFile(delete=False)
            f.write('\n'.join(nl))
            f.close()
            args = ['python',
                    self.rwr_script,
                    network_file,
                    f.name,
                    '-d', delim]
            output = subprocess.check_output(args)
            results.append([tuple(l.strip().split('\t'))
                            for l in output.strip().split('\n')])
            # os.remove(f.name)
        return results


    def _get_nodes(self, network_file):
        nodes = set()
        with open(network_file, 'r') as f:
            for line in f:
                l = line.strip().split()
                assert len(l) == 3
                n1, n2, w = l
                nodes.update((n1, n2))
        return nodes

