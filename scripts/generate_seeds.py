"""
Generate seeds from an edge list
"""
import sys

def main(argv):

    try:
        fp = open(argv[1], 'r')
    except IOError:
        sys.exit('Could not open file: {}'.format(argv[1]))

    output_dir = argv[2]
    nodelist_file = argv[3]

    node_list = set()

    for line in fp.readlines():
        node_list.add(line.rstrip().split()[0])
        node_list.add(line.rstrip().split()[1])

    node_list = list(node_list)

    nodelist_fp = open(nodelist_file, 'w')
    for idx, node in enumerate(node_list):
        nodelist_fp.write(node + '\n')
        output_filename = '{}/seed_{}.txt'.format(output_dir, idx)
        output_fp = open(output_filename, 'w')
        output_fp.write(node)
        output_fp.close()

    nodelist_fp.close()
    fp.close()


if __name__ == '__main__':
    main(sys.argv)
