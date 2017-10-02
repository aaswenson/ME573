import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from problem1 import make_b_array
from solvers import gaussian_elimination, thomas_solver

def make_A_b(n):
     e = np.ones((n, 1))
     data = np.hstack((e, -2*e, e)).T
     diags = np.array((-1, 0, 1))

     sp_array = scipy.sparse.spdiags(data, diags, n, n).toarray()
     
     B = make_b_array(n)

     return sp_array, B

def scaling_analysis():
    gauss_ts = []
    thoms_ts = []
    n = []
    for power in range(2,8):
        print power
        n.append(2 ** power)
        sp_array, B = make_A_b(2 ** power)
        gauss_ts.append(gaussian_elimination(sp_array, B)[2])
        thoms_ts.append(thomas_solver(sp_array, B)[1])
    
    plt.plot(n, gauss_ts)
    plt.plot(n, thoms_ts)
    plt.show()


if __name__ == '__main__':
    
    scaling_analysis()
