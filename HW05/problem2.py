import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from solvers import solve_matrix, gaussian_elimination, thomas_solver
import problem1 as p1

# some constants
k = 1e-3
dx = 0.05
default_dt = 1
t_end = 1e2
domain = (-3, 3)

def make_A_back(n, dx, dt):
    """Produce sparse diagonal matrix for algorithm testing.
    """
    a1 = (1 / dt) + (2 * k / (dx ** 2))
    b1 = - k / (dx ** 2)
    c1 = b1
    
    e = np.ones((n, 1))
    data = np.hstack((c1*e, a1*e, c1*e)).T
    diags = np.array((-1, 0, 1))
    sp_array = scipy.sparse.spdiags(data, diags, n, n).toarray()

    return sp_array

def make_b_array_back(n, dt):
    """ Produce b value array based on problem intial conditions.

             { U = 2  if |x| < 1
    f(x,0) = { 
             {   0    if |x| > 1
    """
    array = np.linspace(-3, 3, n)
    for i, x in enumerate(array[1:-1], start=1):
        if abs(x) < 1:
            array[i] = 2 * (1 / dt)
        else:
            array[i] = 0
    array[0] = 0
    array[n-1] = 0

    return array

def run_btcs_calc(dx, dt):
    
    f = []
    times = np.linspace(0, t_end, t_end / dt)
    n = int( (domain[1] - domain[0]) / dx)
    x = np.linspace( domain[0], domain[1], n)
    b_init = make_b_array_back(n, dt)
    A = make_A_back(n, dx, dt)

    f.append(solve_matrix(A, b_init)[0])

    for i in range(0, len(times)):

        new_b = np.concatenate(([f[i][0]], f[i][1:-1] * (1/dt) ,\
                [f[i][len(f[i])-1]]))
        f.append(solve_matrix(A, new_b)[0])

    return f, x

def get_analytic_solution(x, t):
    fa = np.zeros(len(x))
    for i, x_val in enumerate(x):
        fa[i] = ( 1 ) * ( math.erf( (1-x_val) / (2 * math.sqrt(k * t))) - \
                math.erf(-(1+x_val) / (2 * math.sqrt(k * t) )) )

    return fa

def run_comparison(dt=default_dt):
    time_idx = int(t_end / dt)
    # run ftcs calcs
    f_forward, x = p1.run_ftcs_calc(dx, dt)
    # run btcs calcs
    f_back, x = run_btcs_calc(dx, dt)
    
    # get analytical solution
    fa = get_analytic_solution(x, 1e2)
    
    analytic = (x, fa)
    forward = (x, f_forward[time_idx])
    back = (x, f_back[time_idx])

    return [analytic, forward, back]

def plot(data):

    plt.plot(data[0][0], data[0][1], label='Analytic')
    plt.plot(data[1][0], data[1][1], label='FTCS')
    plt.plot(data[2][0], data[2][1], label='BTCS')
    plt.legend()
    plt.title("Backward in Time, Central in Space Scheme")
    plt.xlabel('x [-]')
    plt.ylabel('f [-]')
    plt.savefig("problem2.png")
    plt.show()

if __name__ == '__main__':
    
    data = run_comparison()
    plot(data)
