import numpy as np
import math
import matplotlib.pyplot as plt

kappa = 0.001
L = 1
N_max = 10000
n_x_steps = 100.0
dx = L / n_x_steps


def eval_f(x):
    """ Evaluate f as a function of x.
    """
    f = 0
    for k in range(1, N_max):
        f += (4.0 / math.pi) * (1.0 / k) * (math.sin(k*math.pi / 4.0) ** 2) *\
        math.sin(2.0 * k * x)
    return x

def eval_u_prime(x):
    """ Evaluate derivative of u as a function of x.
    """
    u_prime = 0
    for k in range(1, N_max):
        u_prime += (math.sin(math.pi * k / 4.0) ** 2) * math.cos(2.0 * k * x)\
        / (2.0 * k)
    return u_prime

def eval_u(x):
    """ Evaluate u as a function of x.
    """
    u = 0
    for k in range(1, N_max):
        u += (-1.0 / math.pi) * (math.sin(math.pi * k / 4.0) ** 2) *\
                math.sin(2.0 * k * x)
    return u

def get_func_data():
    """ Gather data for plotting
    """
    x = []
    f = []
    u = []
    u_prime = []

    for i in range(0, int(n_x_steps)):
        x.append(i * dx)
    for i, val in enumerate(x):
        f.append(eval_f(val))
        u.append(eval_u(val))
        u_prime.append(eval_u_prime(val))
    
    return {'f(x)' : f, 'u(x)' : u, "u'(x)" : u_prime}, x
    
def plot_f(data, x):
    
    for run in data:
        plt.plot(x, data[run], label=run)
    plt.legend()
    plt.title("u(x) estimated with 10,000 Fourier modes")
    plt.xlabel("x [-]")
    plt.ylabel("f(x), u'(x), or u(x) [-]")
    plt.show()

if __name__ == '__main__':

    data, x = get_func_data()
    plot_f(data, x)
