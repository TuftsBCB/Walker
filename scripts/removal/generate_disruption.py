#!/usr/bin/env python
import sys
import argparse
import shared_functions as sh

TOP = 25
CHANGE_THRESHOLD = 20

def generate_summary_mapping(gene_set, tissue_id, unaltered_rank,
                             altered_rank, remove_list):
    #summary_mapping: tissue + removed gene -> (drop count, max drop)
    summary_mapping = []
    for removed_gene in remove_list:
        single_mapping = []
        for gene in gene_set:
            # positive means the rank after the node is removed is greater
            # (i.e. lower in the rank list), negative means rank w/o removing
            # is greater (i.e. it is higher in the rank list)
            try:
                gene_diff = (altered_rank[removed_gene][tissue_id][gene] -
                             unaltered_rank[removed_gene][tissue_id][gene])
                single_mapping.append((gene, gene_diff))
            except KeyError:
                continue
        drop_list = [g for (g, d) in single_mapping
                       if abs(d) > CHANGE_THRESHOLD]
        drop_count = len(drop_list)
        (max_id, max_drop) = sorted(single_mapping, key=lambda x: x[1],
                                    reverse=True)[0]
        summary_mapping.append((removed_gene, drop_count, max_drop, max_id))

    return {rg: (dc, md, id) for (rg, dc, md, id) in summary_mapping}


def write_summary_to_file(mapping, tissue_id, seed_number, id_symbol_map,
                          remove_list):

    list_file = "./seed{}.tis{}.summary".format(seed_number, tissue_id)
    summary_fp = open(list_file, "w")

    # write header
    summary_fp.write("IDsRemoved\tSymbolsRemoved\tDropNumber\tMaxDrop\tMaxSymbol\n")

    for removed_gene in remove_list:
        id_string = ''
        symbol_string = ''
        if isinstance(removed_gene, tuple):
            id_string = ', '.join(removed_gene)
            symbol_string = ', '.join([id_symbol_map[g] for g in removed_gene])
        else:
            id_string = removed_gene
            symbol_string = id_symbol_map[removed_gene]
        id_string = '({})'.format(id_string)
        symbol_string = '({})'.format(symbol_string)

        drop_number = mapping[removed_gene][0]
        max_drop = mapping[removed_gene][1]
        try:
            max_symbol = id_symbol_map[str(mapping[removed_gene][2])]
        except KeyError:
            max_symbol = 'N/A'
        summary_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(id_string, symbol_string,
                                                drop_number, max_drop, max_symbol))

    summary_fp.close()



def main():

    # set up arg parser
    parser = argparse.ArgumentParser()
    parser.add_argument("original_results", help="Original RWR results file")
    parser.add_argument("altered_results",
                        help="RWR results file with nodes removed")
    parser.add_argument("seed_dir",
                        help="Directory where seed files are located")
    parser.add_argument("mapping_file",
                        help="Location of Entrez ID -> protein name mapping file")
    parser.add_argument("-r", "--remove", nargs='+',
                        help="<Optional> Nodes removed in the query experiments")
    parser.add_argument("-rf", "--remove_file", nargs='+',
                        help="<Optional> List of nodes to remove")
    opts = parser.parse_args()

    original_results_dir = opts.original_results
    altered_results_dir = opts.altered_results
    seed_dir = opts.seed_dir
    mapping = opts.mapping_file

    if opts.remove_file:
        remove_list = sh.parse_removal_list_from_file(opts.remove_file[0])
    else:
        remove_list = opts.remove if opts.remove else []

    tissue_list = (14, 20, 22, 33, 46, 51)
    seed_hash = {}
    id_symbol_map = sh.read_mapping(mapping)

    for i in range(1, 5):
        seed_file = "{}/seed{}.asso".format(seed_dir, i)
        seed_hash = sh.read_seed(seed_file, i, seed_hash)

    # generate and print a diff file for each seed (2, 3, 4)
    for seed_number in range(2, 5):
        # we only consider genes in the top 20 in the original network
        # for any one of the tissues that we consider, stored in gene_set
        gene_set = set()
        # gene_rank is passed to read_rwr
        gene_rank = {}
        # altered_rank[seed_number][tissue_id] maps gene numbers to ranks
        # in the tissue being considered
        altered_rank = {}
        unaltered_rank = {}

        for removed_gene in remove_list:
            unaltered_rank[removed_gene] = {}
            altered_rank[removed_gene] = {}
            unaltered_genes = {}
            altered_genes = {}

            for tissue_id in tissue_list:
                prediction_file = "{}/seed.{}.tis.{}.rwr".format(
                        original_results_dir, seed_number, tissue_id)

                if isinstance(removed_gene, tuple):
                    removed_str = '.'.join(removed_gene)
                    altered_pred_file = "{}/seed.{}.tis.{}.remove.{}.rwr".format(
                            altered_results_dir, seed_number, tissue_id, removed_str)
                else:
                    altered_pred_file = "{}/seed.{}.tis.{}.remove.{}.rwr".format(
                            altered_results_dir, seed_number, tissue_id, removed_gene)

                unaltered_genes = sh.generate_sorted_gene_list(prediction_file,
                                                               gene_rank,
                                                               seed_number,
                                                               tissue_id,
                                                               seed_hash)

                altered_genes = sh.generate_sorted_gene_list(altered_pred_file,
                                                             gene_rank,
                                                             seed_number,
                                                             tissue_id,
                                                             seed_hash)

                unaltered_rank[removed_gene][tissue_id] = {g: r
                                        for (g, r) in unaltered_genes}
                altered_rank[removed_gene][tissue_id] = {g: r
                                        for (g, r) in altered_genes}

                # trim to include only the top TOP genes in the diff list
                top_unaltered_genes = zip(*unaltered_genes[:TOP])[0]
                gene_set = gene_set.union(set(top_unaltered_genes))


        # now generate list summarizing each removal
        for tissue_id in tissue_list:
            sorted_mapping = generate_summary_mapping(gene_set, tissue_id,
                                                      unaltered_rank, altered_rank,
                                                      remove_list)

            write_summary_to_file(sorted_mapping, tissue_id, seed_number,
                                   id_symbol_map, remove_list)

if __name__ == "__main__":
    main()
