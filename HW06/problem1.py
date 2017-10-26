import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from solvers import solve_matrix, gaussian_elimination, thomas_solver
from ftcs import make_A_FTCS, make_B_FTCS

# some constants
k = 1e-3
dx = 0.05
dt = 1.5
alpha = k * dt / dx ** 2
t_end = 9
domain = (-3, 3)
n = int((domain[1] - domain[0]) / dx)
x = np.linspace(domain[0], domain[1], n) 
times = np.linspace(0, t_end, int(t_end / dt))

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


def plot(fig, data, savefig):
    plt.figure(fig)
    plt.plot(data[0], data[1], label='Analytic')
    plt.plot(data[0], data[2], label='CN')
    plt.plot(data[0], data[3], label='FTCS')
    plt.legend()
    plt.title("CN vs. FTCS dx = 0.05 dt = 2.0")
    plt.xlabel('x [-]')
    plt.ylabel('f [-]')
#    plt.savefig("{0}.png".format(savefig))
    plt.show()

if __name__ == '__main__':
    
    fa = get_analytic_solution(x, t_end)
    # Set Up CN Stuff
    A = make_A_CN()
    f = np.piecewise(x, [abs(x) < 1], [2])
    B = make_B_CN(f)
    f = solve_matrix(A, B)[0]
    # Set up FTCS stuff
    A_FTCS = make_A_FTCS()
    f_FTCS = np.piecewise(x, [abs(x) < 1], [2])
    B_FTCS = make_B_FTCS(f_FTCS)
    f_FTCS = solve_matrix(A_FTCS, B_FTCS)[0]

    for i in range(0, len(times)):
        # Sovle CN
        B = make_B_CN(f)
        f = solve_matrix(A, B)[0]
        # Solve FTCS
        f_FTCS = solve_matrix(A_FTCS, B_FTCS)[0]
        B_FTCS = make_B_FTCS(f_FTCS)

    data = [x, fa, f, f_FTCS]
    plot(2, data, 'problem1')
