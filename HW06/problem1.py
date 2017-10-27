import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from solvers import solve_matrix, gaussian_elimination, thomas_solver
import ftcs
from constants import *

def make_A_CN():
    """Produce sparse diagonal matrix for algorithm testing.
    """
    e = np.ones((n, 1))
    data = np.hstack((-alpha*e, 2*(alpha+1)*e, -alpha*e)).T
    diags = np.array((-1, 0, 1))
    sp_array = scipy.sparse.spdiags(data, diags, n, n).toarray()
    sp_array[0][0] = 1; sp_array[0][1:] = 0
    sp_array[-1][-1] = 1; sp_array[-1][:-1] = 0

    return sp_array

def make_B_CN(f):
    
    b = np.zeros(n)
    for i, val in enumerate(f[1:-1]):
        b[i] = alpha*f[i-1] + 2*(1-alpha)*f[i] + alpha*f[i+1]
    return b

def get_analytic_solution(x, t):
    fa = np.zeros(n)
    for i, x_val in enumerate(x):
        fa[i] = ( 1 ) * ( math.erf( (1-x_val) / (2 * math.sqrt(k * t))) - \
                math.erf(-(1+x_val) / (2 * math.sqrt(k * t) )) )

    return fa

def get_L_norm(exact, u):
    
    L_norm = max(abs(exact - u))

    return L_norm


def plot(fig, data, savefig):
    plt.figure(fig)
    plt.plot(data[0], data[1], label='Analytic')
    plt.plot(data[0], data[2], label='CN')
    plt.plot(data[0], data[3], label='FTCS')
    plt.legend()
    plt.ylim(0,3)
    plt.title("CN vs. FTCS dx = 0.05 dt = 2.0")
    plt.xlabel('x [-]')
    plt.ylabel('f [-]')
    plt.savefig("./writeup/{0}.png".format(savefig))

if __name__ == '__main__':
    
    fa = get_analytic_solution(x, t_end)
    # Set Up CN Stuff
    A = make_A_CN()
    f = np.piecewise(x, [abs(x) < 1], [2])
    B = make_B_CN(f)
    f = solve_matrix(A, B)[0]


    for i in range(0, len(times)):
        # Sovle CN
        B = make_B_CN(f)
        f = solve_matrix(A, B)[0]
    
    f_FTCS = ftcs.run_ftcs_calc()
    data = [x, fa, f, f_FTCS]
    CN_l_norm = get_L_norm(fa, f)
    FTCS_l_norm = get_L_norm(f_FTCS, f)
    results = "The CN L_inf error is: " + str(CN_l_norm) + "\n" +\
          "The FTCS L_inf error is: " + str(FTCS_l_norm)
    print(results)
    plot(2, data, 'problem1')
