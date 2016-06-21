#!/usr/bin/env python
"""
Script to analyze BRAF removal experiments by generating lists of
top 25 gene differences

"""
import sys
import argparse
import numpy as np
import shared_functions as sh

TOP = 20

def generate_average_mapping(gene_set, tissue_list, unaltered_rank,
                             altered_rank, removed_gene):
    average_mapping = []
    for gene in gene_set:
        diff_list = []
        for tissue_id in tissue_list:
            # positive means the rank after the node is removed is greater
            # (i.e. lower in the rank list), negative means rank w/o removing
            # is greater (i.e. it is higher in the rank list)
            gene_diff = (altered_rank[removed_gene][tissue_id][gene] -
                         unaltered_rank[removed_gene][tissue_id][gene])
            diff_list.append(gene_diff)
        diff_a = np.array(diff_list)
        average_mapping.append((gene, np.mean(diff_a), np.std(diff_a)))

    sorted_mapping = sorted(average_mapping, key=lambda x: abs(x[1]),
                            reverse=True)
    return sorted_mapping


def write_list_to_file(mapping, tissue_id, seed_number, id_symbol_map,
                       removed_gene=None):

    if not removed_gene:
        list_file = "./seed{}.avglist".format(seed_number)
    else:
        if isinstance(removed_gene, tuple):
            removed_str = '.'.join(removed_gene)
            list_file = "./seed{}.remove.{}.avglist".format(
                        seed_number, removed_str)
        else:
            list_file = "./seed{}.remove.{}.avglist".format(
                        seed_number, removed_gene)
    seed_fp = open(list_file, "w")

    # write header
    seed_fp.write("GeneID\tGeneSymbol\tAvgDiff\tStdev\n")

    for (gene, avg, stdev) in mapping:
        try:
            seed_fp.write("{}\t{}\t{:.3f}\t{:.3f}\n".format(gene,
                id_symbol_map[gene], avg, stdev))
        except KeyError:
            # if it's not in our mapping, just write the id number, we can
            # convert it later if necessary (this shouldn't happen often)
            seed_fp.write("{}\tN\A\t{}\t{}\t{}\n".format(gene, diff,
                                                         abs_rank, rem_rank))

    seed_fp.close()

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

        for removed_gene in remove_list:
            sorted_mapping = generate_average_mapping(gene_set, tissue_list,
                                                      unaltered_rank, altered_rank,
                                                      removed_gene)

            write_list_to_file(sorted_mapping, tissue_id, seed_number,
                               id_symbol_map, removed_gene)

if __name__ == "__main__":
    main()
