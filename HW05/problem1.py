import math
import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt

from solvers import solve_matrix, gaussian_elimination, thomas_solver

# some constants
k = 1e-3
default_dt = 0.01
dx_vals = [0.5, 0.25, 0.05]
t_end = 1e2
domain = (-3, 3)

def make_A(n, alpha):
    """Produce sparse diagonal matrix for algorithm testing.
    """
    e = np.ones((n, 1))
    data = np.hstack((alpha * e, (1 - 2*alpha)*e, alpha * e)).T
    diags = np.array((-1, 0, 1))
    sp_array = scipy.sparse.spdiags(data, diags, n, n).toarray()
    sp_array[0][0] = 1
    sp_array[0][1:] = 0
    sp_array[-1][-1] = 1
    sp_array[-1][:-1] = 0
    return sp_array

def make_b_array(n):
    """ Produce b value array based on problem intial conditions.

             { U = 2  if |x| < 1
    f(x,0) = { 
             {   0    if |x| > 1
    """
    array = np.linspace(-3, 3, n)
    for i, x in enumerate(array[1:-1], start=1):
        if abs(x) < 1:
            array[i] = 2
        else:
            array[i] = 0
    array[0] = 0
    array[n-1] = 0

    return array

def run_ftcs_calc(dx, dt=default_dt):
    
    f = []
    times = np.linspace(0, t_end, t_end / dt)
    alpha = k * dt / dx ** 2
    n = int( (domain[1] - domain[0]) / dx)
    x = np.linspace( domain[0], domain[1], n)
    b_init = make_b_array(n)
    A = make_A(n, alpha)
    # set initial conditions and begin forward march
    f.append(b_init)
    for i in range(1, len(times)):
        f1 = f[i-1].dot(A)
        f.append(f1)
     
    return f, x

def get_analytic_solution(x, t):
    fa = np.zeros(len(x))
    for i, x_val in enumerate(x):
        fa[i] = ( 1 ) * ( math.erf( (1-x_val) / (2 * math.sqrt(k * t))) - \
                math.erf(-(1+x_val) / (2 * math.sqrt(k * t) )) )

    return fa

def run_comparison(dt=default_dt):
    
    x_data = []
    f_data = []
    time_idx = int(t_end / dt)
    for dx in dx_vals:
        f, x = run_ftcs_calc(dx)
        f_data.append(f[-1])
        x_data.append(x)
    fa = get_analytic_solution(x, 1e2)
    x_data.append(x)
    f_data.append(fa)

    return x_data, f_data

def plot(x_data, f_data):

    plt.plot(x_data[0], f_data[0], label='0.5')
    plt.plot(x_data[1], f_data[1], label='0.25')
    plt.plot(x_data[2], f_data[2], label='0.05')
    plt.plot(x_data[3], f_data[3], label='Analytic')
    plt.legend()
    plt.title("Forward in Time, Central in Space Scheme")
    plt.xlabel('x [-]')
    plt.ylabel('f [-]')
    plt.savefig("./write_up/problem1.png")
    plt.show()

if __name__ == '__main__':
    
    x_data, f_data = run_comparison()
    plot(x_data, f_data)
