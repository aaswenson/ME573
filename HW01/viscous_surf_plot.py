'''
Plot the x,y profile of viscous forces for an arbitrary fluid in a square duct.
This module will generate the data and produce the 3D surface plot.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

a = 0.1
b = a

step_size = a / (100 * 1.0)

# Make data.
X = np.arange(-a, a, step_size)
Y = np.arange(-b, b, step_size)
X, Y = np.meshgrid(X, Y)

Z = (-7 / (2 * 1.0))*( (1 - X**2/ (a**2) * 1.0)/(b**2 * 1.0) + (1 - Y**2/(b**2 * 1.0))/(a**2 * 1.0))
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(np.min(Z),np.max(Z))
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
