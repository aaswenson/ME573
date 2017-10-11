import matplotlib.pyplot as plt
import math
import numpy as np

def eval_f_prime(x):
    """Evaluate derivative of given f(x)
    """
    f_prime = (5 * math.cos(5*x) / x ** 3) - (3 * math.sin(5*x) / x ** 4)
    return f_prime

def eval_f(x):
    """Evaluate given f(x)
    """
    fx = math.sin(5*x) / (x ** 3)
    return fx

def get_f_f_prime(domain, dx):
    """Set domain and collect, f(x), f'(x) data.
    """
    n = int((domain[1] - domain[0]) / dx)
    x = np.linspace(domain[0], domain[1], n)
    fx = np.zeros(n)
    f_prime = np.zeros(n)
    for idx, val in enumerate(x):
        fx[idx] = eval_f(val)
        f_prime[idx] = abs(eval_f_prime(val))
    
    return [x, fx, f_prime]

def forward_diff(fx, dx):
    """Approximate f'(x) using forward difference numerical scheme.
    """
    f_prime = []
    for idx in range(0, len(fx) - 1):
        f_prime.append((fx[idx + 1] - fx[idx]) / dx)
    
    return f_prime

def second_order_cent(fx, dx):
    """Approximate f'(x) using 2nd order numerical scheme.
    """
    f_prime = []
    for idx in range(1, len(fx) - 1):
        f_prime.append((fx[idx+1] - fx[idx-1]) / (2 * dx))
    return f_prime

def fourth_order_cent(fx, dx):
    """Approximate f'(x) using 4th order numerical scheme.
    """
    f_prime = []
    for idx in range(2, len(fx) - 2):
        f_prime.append((fx[idx-2] - 8*fx[idx-1] + 8*fx[idx+1] - fx[idx+2]) /\
                (12 * dx))
    return f_prime

def compute_num_error(approx, exact, stencil_idx=(None, None)):
    """Compute numerical error of f'(x) approximation.
    """
    
    approx_array = np.array(approx)
    exact_array = np.array(exact[stencil_idx[0]:stencil_idx[1]])

    error = abs(np.subtract(approx_array, exact_array)) / sum(exact_array)

    return error

def compare_numerical_schemes(domain, dx):
    """Compare the numerical approximations of f'(x)
    """
    x, fx, exact = get_f_f_prime(domain, dx)
    er1 = np.zeros(len(exact))
    er2 = np.zeros(len(exact))
    er3 = np.zeros(len(exact))
    # calculate approximations of f'(x)
    first = forward_diff(fx, dx)
    second = second_order_cent(fx, dx)
    fourth = fourth_order_cent(fx, dx)
    # calculate error of approximations
    er1 = compute_num_error(first, exact, (0, -1))
    er2 = compute_num_error(second, exact, (1, -1))
    er4 = compute_num_error(fourth, exact, (2, -2))

    return er1, er2, er4, x

def error_decay(domain):
    
    er1 = []
    er2 = []
    er4 = []
    dx = [0.0001, 0.0005, 0.001, 0.01, 0.05]
    x_eval = 0.2
    f_prime_exact = eval_f_prime(x_eval)
    print f_prime_exact
    for i, h in enumerate(dx):
        # get fx
        x, fx, exact = get_f_f_prime(domain, h)
        # get list index
        idx = int((x_eval - domain[0]) / h) - 1 
        # calculate approximations of f'(x)
        first = forward_diff(fx, h)[idx]
        second = second_order_cent(fx, h)[idx]
        fourth = fourth_order_cent(fx, h)[idx]
        err1 = abs(first - f_prime_exact) / abs(f_prime_exact)
        err2 = abs(second - f_prime_exact) / abs(f_prime_exact)
        err4 = abs(fourth - f_prime_exact) / abs(f_prime_exact)
        print fourth
        print first
        er1.append(err1)
        er2.append(err2)
        er4.append(err4)
    return er1, er2, er4, dx


# plotting functions follow this comment

def plot_parta(data, title):
    """Plot f(x), f'(x) from part A
    """
    plt.figure(1)
    plt.semilogy(data[0], data[1], label='f(x)')
    plt.semilogy(data[0], data[2], label="f'(x)")
    plt.legend()
    plt.xlabel('x [-]')
    plt.ylabel("log(f(x), f'(x)) [-]")
    plt.title(title)
    plt.savefig("prob4_partA.png")
    plt.show()

def plot_partb(dx, data, title):
    plt.figure(2)
    plt.plot(dx[:-1], data[0], label='1st Order')
    plt.plot(dx[1:-1], data[1], label='2nd Order')
    plt.plot(dx[2:-2], data[2], label='4th Order')
    plt.legend()
    plt.xlabel('x [-]')
    plt.ylabel("f'(x)) [-]")
    plt.title(title)
    plt.savefig("prob4_partB.png")
    plt.show()

def plot_partc(h, data, title):
    plt.figure(3)
    plt.loglog(h, data[0], label='1st Order')
    plt.loglog(h, data[1], label='2nd Order')
    plt.loglog(h, data[2], label='4th Order')
    plt.legend()
    plt.xlabel('x [-]')
    plt.ylabel("log(error [-]")
    plt.title(title)
    plt.savefig("prob4_partC.png")
    plt.show()
    # plot error linear
    plt.figure(4)
    plt.plot(h, data[0], label='1st Order')
    plt.plot(h, data[1], label='2nd Order')
    plt.plot(h, data[2], label='4th Order')
    plt.legend()
    plt.xlabel('x [-]')
    plt.ylabel("log(error [-]")
    plt.title(title + 'Linear Plot')
    plt.savefig("prob4_partC_lin.png")
    plt.show()

if __name__ == '__main__':
    
    domain = (0.1, 0.4)
    dx = 0.005
    data = get_f_f_prime(domain, dx)
    plot_parta(data, 'Part A')
    e1, e2, e3, x = compare_numerical_schemes(domain, dx)
    plot_partb(x, [e1, e2, e3], 'Part B')
    er1, er2, er4, h_vals = error_decay(domain)
    plot_partc(h_vals, [er1, er2, er4], 'Part C')
