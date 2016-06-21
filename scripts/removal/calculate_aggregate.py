"""
Script to calculate aggregate stats for dual removal results files

"""
import sys
import argparse
import numpy as np
import shared_functions as sh

def calculate_statistics(input_file):
    diff_list = []
    try:
        input_fp = open(input_file, 'r')
    except IOError:
        sys.exit("Error opening file {}".format(input_file))

    # get rid of headers
    input_fp.readline()

    # get stats from each subsequent line
    for idx, line in enumerate(input_fp.readlines()):
        diff = line.strip().split('\t')[2]
        diff_list.append(int(diff))
        if idx == 0:
            input_max = int(diff)

    diff_a = np.array(diff_list)

    return {
        'Max' : input_max,
        'Mean': np.mean(diff_a),
        'Median': np.median(diff_a),
        'Stdev': np.std(diff_a)
    }


def write_stats_to_file(seed_number, stats, tissue_list, seed_list):
    # write top 5 highest (i.e. most negative):
    #  - max values
    #  - mean values
    #  - median values
    # and list each with stdev within the file

    output_file = "./seed{}.liststats".format(seed_number)
    output_fp = open(output_file, "w")


    max_list = []
    mean_list = []
    median_list = []
    stdev_list = []
    for tissue in tissue_list:
        for gene in seed_list:
            filename = "seed{}.tis{}.remove.{}.{}.toplist".format(
                            seed_number, tissue, 673, gene)
            max_list.append((filename, stats[tissue][gene]['Max']))
            mean_list.append((filename, stats[tissue][gene]['Mean']))
            median_list.append((filename, stats[tissue][gene]['Median']))
            stdev_list.append((filename, stats[tissue][gene]['Stdev']))

    output_fp.write("Top 10 max (absolute) differences:\n")
    for filename, max_value in sorted(max_list, key=lambda x: x[1])[:10]:
        output_fp.write("{}\t{}\n".format(filename, max_value))
    output_fp.write("\n")

    output_fp.write("Top 10 mean (absolute) differences:\n")
    for filename, mean_value in sorted(mean_list, key=lambda x: x[1])[:10]:
        output_fp.write("{}\t{}\n".format(filename, mean_value))
    output_fp.write("\n")

    output_fp.write("Top 10 median (absolute) differences:\n")
    for filename, median_value in sorted(median_list, key=lambda x: x[1])[:10]:
        output_fp.write("{}\t{}\n".format(filename, median_value))
    output_fp.write("\n")

    output_fp.write("Top 10 smallest stdev (absolute) differences:\n")
    for filename, stdev_value in sorted(stdev_list, key=lambda x: x[1])[:10]:
        output_fp.write("{}\t{:.3f}\n".format(filename, stdev_value))


def main():
    # set up arg parser
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir",
                        help="Directory where results files are located")
    parser.add_argument("seed_dir",
                        help="Directory where seed files are located")
    parser.add_argument("mapping_file",
                        help="Location of Entrez ID -> protein name mapping file")
    opts = parser.parse_args()


    results_dir = opts.results_dir
    seed_dir = opts.seed_dir
    mapping = opts.mapping_file

    tissue_list = (14, 20, 22, 33, 46, 51)
    seed2_list = (369, 5894, 5604, 5605, 1326, 5594)
    seed_hash = {}
    stats_hash = {}
    for i in xrange(1, 5):
        seed_file = "{}/seed{}.asso".format(opts.seed_dir, i)
        seed_hash = sh.read_seed(seed_file, i, seed_hash)

    # calculating stats per seed for now
    for seed_number in xrange(2, 5):
        stats_hash[seed_number] = {}

        for tissue_id in tissue_list:
            stats_hash[seed_number][tissue_id] = {}
            for removed_gene in seed2_list:
                stats_hash[seed_number][tissue_id][removed_gene] = {}
                input_file = "{}/seed{}.tis{}.remove.{}.{}.toplist".format(
                      opts.results_dir, seed_number, tissue_id, 673, removed_gene)
                stats_hash[seed_number][tissue_id][removed_gene] = calculate_statistics(input_file)
        write_stats_to_file(seed_number, stats_hash[seed_number], tissue_list, seed2_list)


if __name__ == "__main__":
    main()
