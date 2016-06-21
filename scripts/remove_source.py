"""
Script for removing source node data from results with restart

Copies each file (less the source node) to a new directory, specified
in command line args.

"""

import os
import sys

def main():

    if len(sys.argv) != 3:
        sys.exit("Usage: python {} <source_dir> <dest_dir>".format(sys.argv[0]))

    gene_files = [f for f in os.listdir(sys.argv[1]) if f[:4] == 'gene']

    for filename in gene_files:
        source_id = filename.split('.')[1]

        source_filename = "{}/{}".format(sys.argv[1], filename)
        dest_filename = "{}/{}".format(sys.argv[2], filename)
        in_f = open(source_filename, 'r')
        out_f = open(dest_filename, 'w')

        for line in in_f:
            if line.split()[0] != source_id:
                out_f.write(line)

if __name__ == '__main__':
    main()
