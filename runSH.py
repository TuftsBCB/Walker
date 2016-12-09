#!/usr/bin/env python
import sys
import os

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

ppi = sys.argv[1]
seed_dir = sys.argv[2]
nodelist_file = sys.argv[3]
node_list = get_node_list(sys.argv[3])
output_script = 'run_rwr'

# open script file to write

script_no = 0
script_fp = ''
for idx, gene in enumerate(node_list):

    # for each step, write a command
    if idx % 1000 == 0:
        if script_fp: script_fp.close()
        try:
            script_fp = open(output_script + '{}.sh'.format(script_no), "w")
            script_no += 1
        except IOError:
            sys.exit("Error opening file {}".format(output_script))

    seed_file = '{}/seed_{}.txt'.format(seed_dir, idx)
    output_file = "/r/bcb/TissueSpecificBRAF/string_results/seed.{}.rwr".format(idx)
    command = "python matrix_main.py {} {} -n {} > {}\n".format(
            ppi, seed_file, nodelist_file, output_file)
    script_fp.write(command)

script_fp.close()

