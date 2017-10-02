import numpy as np
import scipy.linalg
import math

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

def solve_matrix(a, b):
    """Solve Ax \ b using numpy's built-in matrix solver.
    """

    x = np.linalg.lstsq(a, b)
    res = np.linalg.norm(np.subtract(a.dot(x[0]), b))
    
    return x[0], res

def gaussian_elimination(A, B, n):
    """Solve the matrix via Gaussian elimination algorithm.
    """
    # LU decomposition with pivot
    pl, u = scipy.linalg.lu(A, permute_l=True)
    # forward substitution to solve for Ly = B
    y = np.zeros(B.size)
    for m, b in enumerate(B.flatten()):
        y[m] = b
        # skip for loop if m == 0
        if m:
            for n in xrange(m):
                y[m] -= y[n] * pl[m,n]
        y[m] /= pl[m, m]
 
    x = np.zeros(B.size)
    lastidx = B.size - 1  # last index
    x[0] = 0
    for midx in xrange(B.size):
        m = B.size - 1 - midx  # backwards index
        x[m] = y[m]
        if midx:
            for nidx in xrange(midx):
                n = B.size - 1  - nidx
                x[m] -= x[n] * u[m,n]
        x[m] /= u[m, m]
        
    res = np.linalg.norm(np.subtract(A.dot(x), B))
    return x, res

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
    g_x, g_res = gaussian_elimination(coeffs, b_array, n)
    
    L_1 = compute_diff(x, g_x)
    results = "The residual of the matrix solver is {0}.\n\
The residual of the Gaussian elimination algorithm is {1}.\n\
The difference between the two algorithms is {2}."\
               .format(res, g_res, L_1)
    print(results)
