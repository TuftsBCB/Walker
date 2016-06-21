import sys

def read_mapping(mapping):
    """ Read (gene ID) -> (gene symbol) mapping file into a dict """
    mapping_hash = {}
    try:
        fp = open(mapping, "r")
    except IOError:
        sys.exit("Error opening file {}".format(mapping))

    for line in fp.readlines():
        info = line.split()
        mapping_hash[info[0]] = info[1]

    fp.close()
    return mapping_hash


def read_seed(seed_file, i, seed_hash):
    """ Read seed file into a list containing sublists representing each seed """
    try:
        fp = open(seed_file, "r")
    except IOError:
        sys.exit("Error opening file {}".format(seed_file))

    seed_hash[i] = []
    for line in fp.readlines():
        info = line.split()
        seed_hash[i].append(info[1])

    fp.close()
    return seed_hash


def read_top_lists(seed_number, seed_hash, pred_dir):
    tissue_list = (14, 20, 22, 33, 46, 51)
    top_lists = {}
    for tissue in tissue_list:
        top_list = []
        filename = "{}/seed.{}.tis.{}.rwr".format(
                        pred_dir, seed_number, tissue)

        try:
            fp = open(filename, "r")
        except IOError:
            sys.exit("Error opening file {}".format(filename))

        # get top 20 genes from file, excepting genes in the
        # current seed
        gene_counter = 1
        while gene_counter < 21:
            current_gene = fp.readline().strip()
            if current_gene not in seed_hash[seed_number]:
                top_list.append(current_gene)
                gene_counter += 1

        top_lists[tissue] = top_list

    return top_lists


def write_top_lists(top_lists, seed_number, id_symbol_map):
    diff_file = "./seed{}.toprank".format(seed_number)
    tissue_list = (14, 20, 22, 33, 46, 51)
    list_fp = open(diff_file, "w")

    # write the header
    list_fp.write("Seed {} Ranking\n\n".format(seed_number))
    for tissue_id in tissue_list:
        list_fp.write("\t{}".format(tissue_id))
    list_fp.write("\n")

    # write rank lists
    for rank in xrange(20):
        for tissue in tissue_list:
            current_gene_id = top_lists[tissue][rank]
            list_fp.write("{}\t\t".format(id_symbol_map[current_gene_id]))
        list_fp.write("\n")


def main(argv):

    if len(sys.argv) != 4:
        sys.exit("python {} <prediction dir> <mapping file>".format(sys.argv[0]))

    pred_dir = argv[1]
    seed_dir = argv[2]
    seed_hash = {}
    mapping = argv[3]
    id_symbol_map = read_mapping(mapping)

    for i in xrange(1, 5):
        seed_file = "{}/seed{}.asso".format(seed_dir, i)
        seed_hash = read_seed(seed_file, i, seed_hash)

    for seed_id in xrange(1, 5):
        top_lists = read_top_lists(seed_id, seed_hash, pred_dir)
        write_top_lists(top_lists, seed_id, id_symbol_map)


if __name__ == "__main__":
    main(sys.argv)
