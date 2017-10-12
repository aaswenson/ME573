import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from solvers import solve_matrix, gaussian_elimination, thomas_solver

# some constants
k = 1e-3
dx = [0.5, 0.25, 0.05]
dt = 0.01

def make_A_b(n, alpha):
    """Produce sparse diagonal matrix for algorithm testing.
    """
    e = np.ones((n, 1))
    data = np.hstack((-alpha * e, (1 + 2*alpha)*e, -alpha * e)).T
    diags = np.array((-1, 0, 1))
    sp_array = scipy.sparse.spdiags(data, diags, n, n).toarray()

    return sp_array

def make_b_array(n):
    """ Produce b value array based on problem intial conditions.

             { U = 2  if |x| < 1
    f(x,0) = { 
             {   0    if |x| > 1
    """
    array = np.linspace(-3, 3, n)
    for i, x in enumerate(array):
        if abs(x) < 1:
            array[i] = 2
        else:
            array[i] = 0
    
    return array

if __name__ == '__main__':

    
    test = make_b_array(10)
    testA = make_A_b(10, 2)

    print test 
    print testA

    f = solve_matrix(testA, test)[0]

    print f
