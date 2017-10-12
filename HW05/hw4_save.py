import numpy as np
import scipy.linalg
import math
from solvers import gaussian_elimination, solve_matrix

def minij(n):
    """Produce a min(i, j) array
    """
    array = np.zeros((n, n))
    for x_idx in range(0, n):
        for y_idx in range(0, n):
            
            array[x_idx][y_idx] = min(x_idx, y_idx) + 1

    return array

def make_b_array(n):
    """ Produce b value array based on 1 + sin(pi * x)
    """

    array = np.linspace(0, 1, n)

    for i, x in enumerate(array):
        array[i] = 1 + math.sin(math.pi * x)

    return array



def compute_diff(mat_x, g_x):
    """Compute difference between Mat solver (Ax \ b) and the Gaussian
    Elimination algorithm
    """
    diff = np.subtract(mat_x, g_x)
    L_1 = abs(sum(diff)) / mat_x.size

    return L_1

if __name__ == '__main__':
    n = 500
    coeffs = minij(n)
    b_array = make_b_array(n)
    x, res = solve_matrix(coeffs, b_array)
    g_x, g_res, t= gaussian_elimination(coeffs, b_array)
    
    L_1 = compute_diff(x, g_x)
    results = "The residual of the matrix solver is {0}.\n\
The residual of the Gaussian elimination algorithm is {1}.\n\
The difference between the two algorithms is {2}."\
               .format(res, g_res, L_1)
    print(results)
