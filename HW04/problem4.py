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

def forward_diff(fx, dx, n):
    """Approximate f'(x) using forward difference numerical scheme.
    """
    f_prime = np.zeros(fx.size)
    for idx in range(0, fx.size - 1):
        f_prime[idx] = (fx[idx + 1] - fx[idx]) / dx
    return f_prime

def second_order_cent(fx, dx, n):
    """Approximate f'(x) using 2nd order numerical scheme.
    """
    f_prime = np.zeros(fx.size)
    for idx in range(1, fx.size - 1):
        f_prime[idx] = (fx[idx+1] - fx[idx-1]) / (2 * dx)
    return f_prime

def fourth_order_cent(fx, dx, n):
    """Approximate f'(x) using 4th order numerical scheme.
    """
    f_prime = np.zeros(fx.size)
    for idx in range(2, fx.size - 2):
        f_prime[idx] = (fx[idx-2] - 8*fx[idx-1] + 8*fx[idx+1] - fx[idx+2]) /\
                (12 * dx)
    return f_prime

def compute_num_error(approx, exact):
    """Compute numerical error of f'(x) approximation.
    """
    
    error = abs(approx - exact)/(sum(exact))

    return error

def compare_numerical_schemes(domain):
    """Compare the numerical approximations of f'(x)
    """
    dx = np.linspace(1e-2, 1e-5)
    for i, dx in enumerate(dx):
        n = (domain[1] - domain[0]) / dx
        exact = get_f_f_prime(domain, dx)[2]
        print exact
        first_order = forward_diff(exact, dx, n)
        second_order = second_order_cent(exact, dx, n)
        fourth_order = fourth_order_cent(exact, dx, n)
        
        first_error = compute_num_error(first_order, exact)
        second_error = compute_num_error(second_order, exact)
        fourth_error = compute_num_error(fourth_order, exact)

    return dx, [first_error, second_error, fourth_error]

# plotting functions follow this comment

def plot_parta(data, title):
    """Plot f(x), f'(x) from part A
    """
    plt.figure()
    plt.semilogy(data[0], data[1], label='f(x)')
    plt.semilogy(data[0], data[2], label="f'(x)")
    plt.legend()
    plt.xlabel('x [-]')
    plt.ylabel("log(f(x), f'(x)) [-]")
    plt.title(title)
#    plt.show()

def plot_partb(dx, data, title):
    plt.figure() 
    plt.plot(dx, data[0])
    plt.plot(dx, data[1])
    plt.plot(dx, data[2])
    plt.show()

if __name__ == '__main__':
    
    domain = (0.1, 0.4)
    dx = 0.0005
    data = get_f_f_prime(domain, dx)
    plot_parta(data, 'Part A')
    step_sizes, results = compare_numerical_schemes(domain)
    plot_partb(step_sizes, results, 'Part B')
