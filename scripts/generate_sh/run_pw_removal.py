#!/usr/bin/env python3
import sys

def read_removal_from_file(remove_file):
    try:
        remove_fp = open(remove_file, 'r')
    except IOError:
        sys.exit('Could not open file: {}'.format(remove_file))

    remove_list = [l.split('\t')[0] for l in remove_fp.readlines()]
    remove_fp.close()
    return remove_list

if len(sys.argv) < 4:
    sys.exit("python {} <ppi> <tissue_prefix> <removal_file> > runAll.sh".format(sys.argv[0]))

ppi = sys.argv[1]
tissue_prefix = sys.argv[2]

prefix_list = tissue_prefix.split("_")
z_threshold = "{}.{}".format(prefix_list[2], prefix_list[3])

python_file = "/r/bcb/TissueSpecificBRAF/python/matrix_main.py"
seed = "/r/bcb/TissueSpecificBRAF/Data/Seed/20150522/seed{}.asso"
tissue_list = (14, 20, 22, 33, 46, 51)
remove_list = read_removal_from_file(sys.argv[3])

for tis in tissue_list:
    low_list = "/r/bcb/TissueSpecificBRAF/Data/Martin/RPKM/T{}/".format(
            z_threshold)
    tis_low = low_list + tissue_prefix + "." + str(tis) + ".dis"

    # open script file to write
    try:
        script_name = "runTis" + str(tis) + ".sh"
        fp = open(script_name, "w")
    except IOError:
        sys.exit("Error opening file {}".format(script_name))

    # for each gene, write a command
    for seed_id in range(2, 5):
        seed_file = seed.format(seed_id)

        for node_to_remove in remove_list:
            output_file = "seed.{}.tis.{}.remove.{}.{}.rwr\n".format(
                        seed_id, tis, 673, node_to_remove)
            command = "python {} {} {} {} -r {} {} > {}".format(
                    python_file, ppi, tis_low, seed_file,
                    673, node_to_remove, output_file)
            fp.write(command)

    fp.close()

    print("sh runTis{}.sh &".format(tis))


