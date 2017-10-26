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
domain = (-2, 2)
n = int((domain[1] - domain[0]) / dx)
x = np.linspace(domain[0], domain[1], n) 
times = np.linspace(0, t_end, int(t_end / dt))

def make_A_back(n, dx, dt):
    """Produce sparse diagonal matrix for algorithm testing.
    """
    e = np.ones((n, 1))
    data = np.hstack((-alpha*e, 2*(alpha+1)*e, -alpha*e)).T
    diags = np.array((-1, 0, 1))
    sp_array = scipy.sparse.spdiags(data, diags, n, n).toarray()
    sp_array[0][0] = 1
    sp_array[0][1:] = 0
    sp_array[-1][-1] = 1
    sp_array[-1][:-1] = 0
    print sp_array 

    return sp_array

def eval_init_f(x):
    """ Produce proper value of f
    """
    if abs(x) < 1:
        return  2 
    else:
        return 0

def make_b_array_back(n, dt):
    """ Produce b value array based on problem intial conditions.

             { U = 2  if |x| < 1
    f(x,0) = { 
             {   0    if |x| > 1
    """
    b_0 = np.zeros(n)
    for i, val in enumerate(x[1:-2], start=1):
        b_0[i] = alpha * eval_init_f(x[i-1])\
                + 2*(1 - alpha) * eval_init_f(x[i])\
                 + alpha * eval_init_f(x[i+1])
    b_0[0] = 0
    b_0[n-1] = 0
    
    return b_0

def update_b_array(n, f_0):
    b_0 = np.linspace(-3, 3, n)
    for i, val in enumerate(f_0[1:-2], start=1):
        b_0[i] = alpha * f_0[i-1]\
                + 2*(1 - alpha) * f_0[i]\
                 + alpha * f_0[i+1]
    b_0[0] = f_0[0]
    b_0[n-1] = f_0[-1]
    
    return b_0

def run_btcs_calc(dx):
    
    f = []
    times = np.linspace(0, t_end, t_end / dt)
    b_init = make_b_array_back(n, dt)
    A = make_A_back(n, dx, dt)
    f.append(solve_matrix(A, b_init)[0])

    for i in range(0, len(times)):
        new_b = update_b_array(n, f[i])
        f.append(solve_matrix(A, new_b)[0])
    
    return f, x

def get_analytic_solution(x, t):
    fa = np.zeros(len(x))
    for i, x_val in enumerate(x):
        fa[i] = ( 1 ) * ( math.erf( (1-x_val) / (2 * math.sqrt(k * t))) - \
                math.erf(-(1+x_val) / (2 * math.sqrt(k * t) )) )

    return fa

def run_comparison():
    time_idx = int(t_end / dt)
    # run btcs calcs
    f_back, x = run_btcs_calc(dx)
    # get analytical solution
    fa = get_analytic_solution(x, 1e2)
    analytic = (x, fa)
    back = (x, f_back[0])

    return [analytic, back]

def plot(data, savefig):

    plt.plot(data[0][0], data[0][1], label='Analytic')
    plt.plot(data[1][0], data[1][1], label='FTCS')
    plt.legend()
    plt.title("FTCS vs. BTCS dx = 0.05 dt = 2.0")
    plt.xlabel('x [-]')
    plt.ylabel('f [-]')
    plt.savefig("{0}.png".format(savefig))
    plt.show()

if __name__ == '__main__':
    
    data1 = run_comparison()
    plot(data1, "problem1")
