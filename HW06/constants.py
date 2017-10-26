import numpy as np

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
