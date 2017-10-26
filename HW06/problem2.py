import math
import numpy as np
import scipy.sparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from solvers import solve_matrix, gaussian_elimination, thomas_solver

# some constants
k = 0.1
dx = 0.05
dy = 0.05
dt = (dx**2 / k) / 4.0 
t_end = 80 * dt
x_domain = (0, 1)
y_domain = (0, 1)
N_x = int((x_domain[1] - x_domain[0]) / dx)
N_y = int((y_domain[1] - y_domain[0]) / dy)
x = np.linspace(x_domain[0], x_domain[1], N_x) 
y = np.linspace(y_domain[0], y_domain[1], N_y) 
times = np.linspace(0, t_end, int(t_end / dt))

alpha = k * dt / dx ** 2

def get_analytic_solution(x, y, t, trunc):

    X, Y = np.meshgrid(x, y)
    Z = 0
    for n in range(1, trunc):
        for m in range(1, trunc):
            Z_num = -120 * ( ((-n)**4 * np.pi**4 * (-1)**n) +\
                             (12 * n**2 * np.pi ** 2 * (-1)**n)\
                             + 24 + (24 * (-1)**(1+n))\
                             *(-2 + (2*(-1)**m) ) )
            Z_num_xy = np.sin(n*X*np.pi)*np.sin(m*Y*np.pi)\
                       * np.exp(-(n**2 + m**2) * np.pi**2 * k * t)
            Z_denom = n**7 * np.pi**10 * m**3

            Z += Z_num * Z_num_xy / Z_denom
    
    return X, Y, Z

if __name__ == '__main__':

    X, Y, Z = get_analytic_solution(x, y, t_end, 100)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

