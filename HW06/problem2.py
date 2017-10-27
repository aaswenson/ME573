import numpy as np
import matplotlib.pyplot as plt

# some important constants 
x_bound = y_bound = 1.
dx = dy = 0.05
k = 0.1
nx, ny = int(x_bound/dx), int(y_bound/dy)
dx2, dy2 = dx*dx, dy*dy
dt = (dx2 / k) / 4.0
t_end = 80 * dt

# set the grid
u0 = np.zeros((nx, ny))
u_exact = np.zeros((nx, ny))
u = np.zeros((nx, ny))

def get_exact(x, y, t, trunc):
    """Get the exact solution at a set t
    """
    Z=0
    for n in range(1, trunc):
        for m in range(1, trunc):
            Z_num = -120 * ( ((-n)**4 * np.pi**4 * (-1)**n) +\
                             (12 * n**2 * np.pi ** 2 * (-1)**n)\
                             + 24 + (24 * (-1)**(1+n))\
                             *(-2 + (2*(-1)**m) ) )
            Z_num_xy = np.sin(n*x*np.pi)*np.sin(m*y*np.pi)\
                       * np.exp(-(n**2 + m**2) * np.pi**2 * k * t)
            Z_denom = n**7 * np.pi**10 * m**3

            Z += Z_num * Z_num_xy / Z_denom
    
    return Z

def get_L_norm(exact, u):
    
    diffs = abs(exact - u)
    l_diffs = []
    for row in diffs:
        l_diffs.append(max(row))
    return max(l_diffs), diffs

# Initial conditions

for i in range(nx):
    for j in range(ny):
        x = i*dx
        y = j*dy
        u0[i,j] = x * (1-x**5) * y * (1-y)
        u_exact[i,j] = get_exact(x, y, t_end, 10)


def do_timestep(u0, u):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + k * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )

    u0 = u.copy()
    return u0, u

u0, u = do_timestep(u0, u)
l_inf_norm, norm_diff_vals = get_L_norm(u_exact, u0)

fig = plt.figure(1)
ax = fig.add_subplot(111)
im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=0, vmax=0.06)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
ax.set_title('2D distribution after 80 time steps using FTCS')
plt.xlabel('x node [-]')
plt.ylabel('y node [-]')
fig.colorbar(im, cax=cbar_ax)
plt.savefig('./writeup/problem2_plot.png')

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.set_title('|f_exact - f_ftcs| Using FTCS')
plt.xlabel('x node [-]')
plt.ylabel('y node [-]')
im = ax.imshow(norm_diff_vals.copy(), cmap=plt.get_cmap('hot'), vmin=0,\
        vmax=l_inf_norm)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
fig.colorbar(im, cax=cbar_ax)
plt.savefig('./writeup/problem2_error.png')

print("The L_infinity error for FTCS is: " + str(l_inf_norm))
