import sys
import numpy as np

def main(argv):
    file_prefix = argv[1]
    num_files = int(argv[2])
    output_filename = argv[3]
    matrix = []

    for idx in range(num_files):
        filename = '{}.{}.rwr'.format(file_prefix, idx)
        try:
            fp = open(filename, 'r')
        except IOError:
            sys.exit('Could not open file: {}'.format(filename))

        matrix.append(np.loadtxt(filename))

    np.savetxt(output_filename, np.array(matrix), fmt='%.10f')

if __name__ == '__main__':
    main(sys.argv)
