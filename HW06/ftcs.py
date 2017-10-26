import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from solvers import solve_matrix, gaussian_elimination, thomas_solver

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

def make_A_FTCS():
    """Produce sparse diagonal matrix for algorithm testing.
    """
    e = np.ones((n, 1))
    data = np.hstack((alpha * e, (1 - 2*alpha)*e, alpha * e)).T
    diags = np.array((-1, 0, 1))
    sp_array = scipy.sparse.spdiags(data, diags, n, n).toarray()
    sp_array[0][0] = 1; sp_array[0][1:] = 0
    sp_array[-1][-1] = 1; sp_array[-1][:-1] = 0

    return sp_array

def make_B_FTCS(f):
    
    b = f
    b[0] = 0
    b[n-1] = 0
     
    return b

def get_analytic_solution(x, t):
    fa = np.zeros(len(x))
    for i, x_val in enumerate(x):
        fa[i] = ( 1 ) * ( math.erf( (1-x_val) / (2 * math.sqrt(k * t))) - \
                math.erf(-(1+x_val) / (2 * math.sqrt(k * t) )) )

    return fa


def plot(data, savefig):

    plt.plot(data[0], data[1], label='Analytic')
    plt.plot(data[0], data[2], label='CN')
    plt.legend()
    plt.title("CN vs. FTCS dx = 0.05 dt = 2.0")
    plt.xlabel('x [-]')
    plt.ylabel('f [-]')
    plt.savefig("{0}.png".format(savefig))
    plt.show()

def run_ftcs():
    A = make_A_FTCS()
    f = np.piecewise(x, [abs(x) < 1], [2])
    B = make_B_FTCS(f)
    f = solve_matrix(A, B)[0]

    for i in range(0, len(times)):
        B = make_B_CN(f)
        f = solve_matrix(A, B)[0]

    return f

if __name__ == '__main__':
    
    A = make_A_FTCS()
    f = np.piecewise(x, [abs(x) < 1], [2])
    B = make_B_FTCS(f)
    f = solve_matrix(A, B)[0]

    for i in range(0, len(times)):
        B = make_B_FTCS(f)
        f = solve_matrix(A, B)[0]

    fa = get_analytic_solution(x, t_end)
    data = [x, fa, f]
    plot(data, 'problem1.png')
