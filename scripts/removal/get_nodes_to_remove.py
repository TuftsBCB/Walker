import sys
import argparse
import shared_functions as sh

TOP = 25

def write_gene_set_to_file(gene_set, id_symbol_map):
    remove_file = "./to_remove.txt"
    remove_fp = open(remove_file, 'w')
    for gene in gene_set:
        remove_fp.write("{}\t{}\n".format(gene, id_symbol_map[gene]))
    remove_fp.close()

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
    remove_list = [673]
    gene_set = set()

    for i in xrange(1, 5):
        seed_file = "{}/seed{}.asso".format(seed_dir, i)
        seed_hash = sh.read_seed(seed_file, i, seed_hash)

    # generate and print a diff file for each seed (2, 3, 4)
    for seed_number in xrange(2, 5):
        # we only consider genes in the top 20 in the original network
        # for any one of the tissues that we consider, stored in gene_set
        # gene_rank is passed to read_rwr
        gene_rank = {}
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
                top_altered_genes = zip(*altered_genes[:TOP])[0]
                gene_set = gene_set.union(set(top_unaltered_genes))
                gene_set = gene_set.union(set(top_altered_genes))

    write_gene_set_to_file(gene_set, id_symbol_map)

if __name__ == "__main__":
    main()
