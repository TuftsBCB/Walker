"""
Version of generateDifferenceTable script for pairwise comparison of genes
between different versions of the graph (i.e. with some nodes removed), within
a single tissue.

"""
import sys
import argparse
import shared_functions as sh

TOP = 100

def generate_sorted_gene_list(prediction_file, gene_rank, seed_no, tissue_id, seed_hash):
    """ Gets a list of (gene, rank) tuples, sorted by rank """
    if seed_no not in gene_rank:
        gene_rank[seed_no] = {}

    gene_rank[seed_no][tissue_id] = sh.read_rwr(prediction_file,
                                                 seed_hash,
                                                 seed_no)
    filtered_genes = [(g, r) for g, r
                      in gene_rank[seed_no][tissue_id].iteritems()]

    sorted_genes = sorted(filtered_genes, key=lambda x: x[1])
    return sorted_genes


def write_list_to_file(diff_list, tissue_id, seed_number, id_symbol_map,
                       removed_gene=None):
    """ Write the calculated diff scores to a toplist file """

    if not removed_gene:
        list_file = "./seed{}.tis{}.toplist".format(seed_number, tissue_id)
    else:
        list_file = "./seed{}.tis{}.remove.{}.toplist".format(
                    seed_number, tissue_id, removed_gene)
    seed_fp = open(list_file, "w")

    # write header
    seed_fp.write("GeneID\tGeneSymbol\tDiff\tNormalRank\tRemoveRank\n")

    for (gene, diff, abs_rank, rem_rank) in diff_list[:25]:
        try:
            seed_fp.write("{}\t{}\t{}\t{}\t{}\n".format(gene, id_symbol_map[gene], diff,
                                                        abs_rank, rem_rank))
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
    opts = parser.parse_args()

    original_results_dir = opts.original_results
    altered_results_dir = opts.altered_results
    seed_dir = opts.seed_dir
    mapping = opts.mapping_file

    tissue_list = (14, 20, 22, 33, 46, 51)
    seed_hash = {}
    id_symbol_map = sh.read_mapping(mapping)

    for i in xrange(1, 5):
        seed_file = "{}/seed{}.asso".format(seed_dir, i)
        seed_hash = sh.read_seed(seed_file, i, seed_hash)

    # generate and print a diff file for each seed (1, 2, 3, 4)
    for seed_number in xrange(2, 5):
        gene_rank = {}

        for tissue_id in tissue_list:
            diff_list = []
            unaltered_genes = {}
            altered_genes = {}

            # have to generate sorted genes for both files
            # then take the top 100, etc.

            # adding loop through seed 2 genes for remove BRAF +
            # seed 2 genes experiment
            seed2_list = (369, 5894, 5604, 5605, 1326, 5594)
            for removed_gene in seed2_list:
                diff_list = []
                prediction_file = "{}/seed.{}.tis.{}.rwr".format(
                        original_results_dir, seed_number, tissue_id)
                altered_pred_file = "{}/seed.{}.tis.{}.remove.{}.rwr".format(
                        altered_results_dir, seed_number, tissue_id,
                        removed_gene)

                unaltered_genes = generate_sorted_gene_list(prediction_file,
                                                            gene_rank,
                                                            seed_number,
                                                            tissue_id,
                                                            seed_hash)

                altered_genes = generate_sorted_gene_list(altered_pred_file,
                                                          gene_rank,
                                                          seed_number,
                                                          tissue_id,
                                                          seed_hash)

                # trim to include only the top TOP genes in the diff list
                top_unaltered_genes = zip(*unaltered_genes[:TOP])[0]
                top_altered_genes = zip(*altered_genes[:TOP])[0]
                unaltered_index = {g: r for (g, r) in unaltered_genes}
                altered_index = {g: r for (g, r) in altered_genes}

                gene_set = set(top_unaltered_genes).union(set(top_altered_genes))

                for gene in gene_set:
                    # positive means the original rank (without removing BRAF)
                    # is greater, negative means rank w/o BRAF is greater
                    gene_diff = unaltered_index[gene] - altered_index[gene]
                    diff_list.append((gene, gene_diff, unaltered_index[gene], altered_index[gene]))

                # sort by absolute value of difference
                diff_list = sorted(diff_list, key=lambda x: abs(x[1]), reverse=True)
                write_list_to_file(diff_list, tissue_id, seed_number,
                                   id_symbol_map, removed_gene)


if __name__ == "__main__":
    main()
