import time
import scipy.linalg
import math
import numpy as np

def solve_matrix(a, b):
    """Solve Ax \ b using numpy's built-in matrix solver.
    """

    x = np.linalg.lstsq(a, b)
    res = np.linalg.norm(np.subtract(a.dot(x[0]), b))
    
    return x[0], res

def gaussian_elimination(A, B):
    """Solve the matrix via Gaussian elimination algorithm.
    """
    start = time.time()
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
    end = time.time()
    res = np.linalg.norm(np.subtract(A.dot(x), B))
    ex_time = (end - start) * 1000 

    return x, res, ex_time

def thomas_solver(A, B):
    start = time.time()
    for k in range(0, B.size - 1):
        i = k + 1
        l = A[i][k] / A[k][k]
        for j in range(k, k + 1):
            A[i][j] -= l * A[k][j]
        B[i] -= l * B[k]
    x = np.zeros(B.size)
    for midx in xrange(B.size):
        m = B.size - 1 - midx  # backwards index
        x[m] = B[m]
        if midx:
            for nidx in xrange(midx):
                n = B.size - 1  - nidx
                x[m] -= x[n] * A[m,n]
        x[m] /= A[m, m]
    end = time.time()
    ex_time = (end - start) * 1000
    return x, ex_time
