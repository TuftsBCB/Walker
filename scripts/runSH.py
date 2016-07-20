#!/usr/bin/env python
import sys

if len(sys.argv) < 5:
    sys.exit("python {} <ppi> <tissuePrefix> <maxStep> <runNumber> > runAll.sh".format(sys.argv[0]))

ppi = sys.argv[1]
# prefix to the low list, i.e. "z_scores_0_1_thres.csv"
tissue_prefix = sys.argv[2]
max_step = int(sys.argv[3])
run_number = sys.argv[4]

prefix_list = tissue_prefix.split("_")
z_threshold = "{}.{}".format(prefix_list[2], prefix_list[3])

python_file = "/r/bcb/TissueSpecificBRAF/python/graph_main.py"
seed = "/r/bcb/TissueSpecificBRAF/Data/Seed/20150522/seed4.asso"
tissue_list = (14, 20, 22, 33, 46, 51)
gene_list = []

try:
    fp = open(seed, "r")
except IOError:
    sys.exit("Error opening file {}".format(seed))

for line in fp.readlines():
    info = line.split()
    gene_list.append(info[1])

fp.close()

for tis in tissue_list:
    low_list = "/r/bcb/TissueSpecificBRAF/Data/Martin/RPKM/T{}/".format(
            z_threshold)
    tis_low = low_list + tissue_prefix + "." + str(tis) + ".dis"

    # open script file to write
    try:
        script_name = "runTis" + str(tis) + ".sh"
        fp = open(script_name, "w")
    except IOError:
        sys.exit("Error opening file {}".format(seed))

    # for each step, write a command
    # for step_i in range(1, max_step + 1):
    for step_i in (3, 5, 8, 10):
        for gene in gene_list:
            output_file = "gene.{}.tis.{}.step.{}.rwr\n".format(
                    gene, tis, step_i)
            command = "python {} {} {} {} {} > {}".format(
                    python_file, ppi, tis_low, gene, step_i, output_file)
            fp.write(command)

    fp.close()

    print("sh runTis{}.sh".format(tis))


