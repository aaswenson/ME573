import numpy as np
import math
import matplotlib.pyplot as plt

kappa = 0.001
L = 1
N_max = 200
n_x_steps = 50.0
dx = L / n_x_steps

def step(x):
    """ Return heaviside eqn value"""

    x_eval = x - (L / 2.0)
    
    return 1 * (x_eval > 0)

def f_x_0(x):

    step_switch = step(x)

    f = x + ((L-x) - x) * step_switch

    return f

def f_x_t(x, t, n):

    fval = 0
    for i in range(1,n + 1):
        C_i = (4 * L / ((math.pi * i) **  2)) * math.sin(i * math.pi * x / L)
        fval += C_i * math.sin(i * math.pi * x / L)*\
                      math.exp(-i ** t)
    return fval

def eval_function(t, N_max, erf=False):

    x_save = []
    f_x_save = []
    for i in range(0, int(n_x_steps)):
        x = i * dx
        if (t == 0) and (erf==False):
            f_x = f_x_0(x)
        else:
            f_x = f_x_t(x, t, N_max)
        x_save.append(x)
        f_x_save.append(f_x)
    
    return x_save, f_x_save

def calc_error(trunc_dat, exact_dat, n_x_steps):

    error = 0

    for i in range(0, int(n_x_steps)):
        err = abs(trunc_dat[i] - exact_dat[i])
        if err > error:
            error = err
    return error


def plot_f(data):
    
    for run in data:
        plt.plot(data[run][0], data[run][1], label=run)
    plt.legend()
    plt.title("f(x,t) Fourier Modes (N-max = 200)")
    plt.xlabel("x [-]")
    plt.ylabel("f(x,t) [-]")
    plt.show()

def plot_err(data):
    
    for run in data:
        plt.semilogy(data[run][0], data[run][1], label=run)
    plt.legend()
    plt.title("Numerical Error is Dependent on Number of Fourier Modes")
    plt.xlabel("N modes [-]")
    plt.ylabel("log(error) [-]")
    plt.show()

if __name__ == '__main__':

    t_coeffs = {'t_0' : 0, 't_med' : 0.1, 't_end' : 0.5}
    data_exact = {'t_0' : [], 't_med' : [], 't_end' : []}
    error  = {'t_0' : [[],[]], 't_med' : [[],[]], 't_end' : [[],[]]}
    
    for time in t_coeffs:
        t = t_coeffs[time]
        x, f_vals = eval_function(t, N_max)
        data_exact[time].append(x)
        data_exact[time].append(f_vals)

        for modes in range(1, 50):
            x, f_vals = eval_function(t, modes, True)
            err = calc_error(f_vals, data_exact[time][1], modes)
            
            error[time][0].append(modes)
            error[time][1].append(err)
    # plot data/error
    plot_f(data_exact)
    plot_err(error)

