import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from solvers import solve_matrix, gaussian_elimination, thomas_solver
from constants import *

def run_ftcs_calc():
    
    # Set Up FTCS matrix approx
    f = np.piecewise(x, [abs(x) < 1], [2])
    new_f = np.zeros(len(f))
    for i in range(0, len(times)):
        for j, val in enumerate(f[1:-1]):
            new_f[j] = alpha * f[j+1] + (1-2*alpha)*f[j] + alpha*f[j-1]
            new_f[0] = 0
            new_f[-1] = 0
        f = new_f
        
    return f
