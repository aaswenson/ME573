import numpy as np
import scipy.linalg
import math

def minij(n):
    """Produce a min(i, j) array
    """
    array = np.zeros((n, n))
    for x_idx in range(0, n):
        for y_idx in range(0, n):
            
            array[x_idx][y_idx] = min(x_idx, y_idx)

    return array

def make_b_array(n):

    array = np.linspace(0, 1, n)

    for i, x in enumerate(array):
        array[i] = 1 + math.sin(math.pi * x)

    return array

def solve_matrix(a, b):

    x = np.linalg.lstsq(a, b)
    res = np.linalg.norm(np.subtract(a.dot(x[0]), b))
    
    return x[0], res

def gaussian_elimination(a, b, n):
    tri = scipy.linalg.lu(a)[2] * a
    print a
    idx = n - 1
    x = np.zeros(n)
    x[idx] = b[idx] / tri[idx][idx]

    for i in reversed(range(0, idx)):
        x[i] = b[i]
        for j in reversed(range(i + 1, idx)):
            x[i] = x[i] - tri[i][j] * x[j]
        x[idx] = x[idx] / tri[idx][idx]

    print x

if __name__ == '__main__':
    n = 5
    coeffs = minij(n)
    b_array = make_b_array(n)
    x_vals, res = solve_matrix(coeffs, b_array)
    gaussian_elimination(coeffs, b_array, n)
    print x_vals
