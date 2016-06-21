"""
Script for making beautiful tables from output

"""

import re
import sys

TISSUE_MAPPING = '/r/bcb/TissueSpecificBRAF/Data/Martin/means.dis.list'

def get_clean_indices(source_p, tissue_hash):
    # first line - map tissue index to tissue name
    indices = source_p.readline().split()
    clean_indices = []
    for index in indices:
        try:
            int(index)
            # remove parentheses and everything in them (this is
            # mostly to prevent headings from getting too long)
            clean_index = re.sub(r"\(.*\)", "", tissue_hash[index])
            clean_indices.append(clean_index)
        except ValueError:
            clean_indices.append(index)

    return clean_indices


def main(argc, argv):

    if argc != 2:
        sys.exit("Usage: python {} <original_table>".format(argv[0]))

    tissue_hash = {}
    with open(TISSUE_MAPPING, 'r') as fp:
        for line in fp.readlines():
            # only split on first tab
            tissue_info = line.split('\t', 1)
            try:
                tissue_hash[tissue_info[0]] = tissue_info[1].strip()
            except IndexError:
                pass

    with open(argv[1], 'r') as source_p:
        dest_file = "{}.clean".format(str(argv[1]))
        with open(dest_file, 'w') as dest_p:
            clean_indices = get_clean_indices(source_p, tissue_hash)
            longest_index = max([len(x) for x in clean_indices])
            format_string = "|".join(['{{:{}}}'.format(longest_index) for
                                       _ in range(len(clean_indices))])
            format_string += '\n'
            dest_p.write(format_string.format(*clean_indices))
            for line in source_p.readlines():
                line_entries = line.split()
                dest_p.write(format_string.format(*line_entries))


if __name__ == '__main__':
    main(len(sys.argv), sys.argv)
