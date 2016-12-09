import sys
import numpy as np

def main(argv):
    matrix_file = argv[1]
    output_filename = argv[2]
    matrix = np.loadtxt(matrix_file)
    def f(x): return 1 / float(x)
    f = np.vectorize(f)
    matrix = f(matrix)

    transformed_matrix = [([0] * len(matrix[0])) for _ in xrange(len(matrix[0]))]
    for i, row in enumerate(matrix):
        for j, col in enumerate(row):
            transformed_matrix[i][j] = matrix[i][j] + matrix[j][i]

    np.savetxt(output_filename, np.array(transformed_matrix), fmt='%.10f')

if __name__ == '__main__':
    main(sys.argv)
