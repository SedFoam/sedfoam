import numpy as np
from pylab import mlab
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import fluidfoam


# Function reading the mesh dimensions
def readXY(sol):
    X, Y, Z = fluidfoam.readmesh(sol)
    return X, Y


# Function returning the bed profile
def depth(sol, t, x, y, xi, yi):
    Nx = np.size(xi,1)
    ybed = np.zeros(Nx)
    if np.mod(t, 1) == 0:
        timename = str(int(t)) + '/'
    else:
        timename = str(t) + '/'
    alpha = fluidfoam.readscalar(sol, timename, 'alpha_a')
    alphai = griddata((x, y), alpha, (xi, yi))
    for j in range(Nx - 1):
        tab = np.where(alphai[:,j+1] > 0.5)
        ybed[j] = yi[np.max(tab),j+1]
    return ybed


# Case information
case_toplot = '../2DPipelineScour/'
time1 = 11
time2 = 18
time3 = 25

# Domain dimensions
D = 0.05
xmin = -0.1
xmax = 0.3
ymin = -0.08
ymax = 0.075

# Number of division for linear interpolation
ngridx = 1500
ngridy = 500

# Interpolation grid dimensions
xinterpmin = -0.1
xinterpmax = 0.35
yinterpmin = -0.09
yinterpmax = 0.015

# Interpolation grid
xi = np.linspace(xinterpmin, xinterpmax, ngridx)
yi = np.linspace(yinterpmin, yinterpmax, ngridy)

# Bed elevation calculation
x1, y1 = readXY(case_toplot)
Xi, Yi = np.meshgrid(xi, yi)

ybed1 = depth(case_toplot, time1, x1, y1, Xi, Yi)
ybed2 = depth(case_toplot, time2, x1, y1, Xi, Yi)
ybed3 = depth(case_toplot, time3, x1, y1, Xi, Yi)

# Experimental data collection
expe_11s = 'DATA/Mao_11s_expe.txt'
expe_18s = 'DATA/Mao_18s_expe.txt'
expe_25s = 'DATA/Mao_25s_expe.txt'

x_expe_11s, y_expe_11s = np.loadtxt(expe_11s, usecols=(0, 1), unpack=True,
                                    delimiter=',')
x_expe_18s, y_expe_18s = np.loadtxt(expe_18s, usecols=(0, 1), unpack=True,
                                    delimiter=',')
x_expe_25s, y_expe_25s = np.loadtxt(expe_25s, usecols=(0, 1), unpack=True,
                                    delimiter=',')

# -------------PLOT SECTION------------- #
# Figure dimensions and parameters
fig_sizex = 27.5
fig_sizey = 37
font_size = 40
figname = 'bed_profiles'
line_style = '-'
line_color = 'C1'
line_width = 4
label_num = 'sedFoam'
marker_size = 16
label_expe = 'Mao [1986]'

# Figure creation
fig = plt.figure(figsize=(fig_sizex, fig_sizey), dpi=100)
plt.rcParams.update({'font.size': font_size})

# Subplots creation
# Subplot 1
plt1 = plt.subplot(3, 1, 1)
circle = plt.Circle((0, 0), radius=0.5, fc='silver', edgecolor='k',
                    zorder=3)
plt.gca().add_patch(circle)
plt1.set_xticklabels([])
plt1.grid()

# Subplot 2
plt2 = plt.subplot(3, 1, 2, sharey=plt1)
circle = plt.Circle((0, 0), radius=0.5, fc='silver', edgecolor='k',
                    zorder=3)
plt.gca().add_patch(circle)
plt2.set_xticklabels([])
plt2.grid()

# Subplot 3
plt3 = plt.subplot(3, 1, 3, sharey=plt1)
circle = plt.Circle((0, 0), radius=0.5, fc='silver', edgecolor='k',
                    zorder=3)
plt.gca().add_patch(circle)
plt3.grid()

plt.subplots_adjust(hspace=0.2)

# Subplot axis
plt1.axis([xmin/D, xmax/D, ymin/D, ymax/D])
plt2.axis([xmin/D, xmax/D, ymin/D, ymax/D])
plt3.axis([xmin/D, xmax/D, ymin/D, ymax/D])
plt1.set_xlabel('')
plt1.set_ylabel('y/D')
plt2.set_ylabel('y/D')
plt3.set_xlabel('x/D')
plt3.set_ylabel('y/D')

# Horizontal line at y = 0
n = np.zeros(2)
nx = np.linspace(xmin/D, xmax/D, 2)
plt1.plot(nx, n-0.5, "k--")
plt2.plot(nx, n-0.5, "k--")
plt3.plot(nx, n-0.5, "k--")

# Numerical results plotting
plt1.plot(xi/D, ybed1/D, label=label_num, linewidth=line_width, ls=line_style,
          color=line_color)
plt2.plot(xi/D, ybed2/D, linewidth=line_width, ls=line_style, color=line_color)
plt3.plot(xi/D, ybed3/D, linewidth=line_width, ls=line_style, color=line_color)

# Experimental data plotting
plt1.plot(x_expe_11s/D, y_expe_11s/D, "ro", label=label_expe,
          markersize=marker_size)
plt2.plot(x_expe_18s/D, y_expe_18s/D, "ro", markersize=marker_size)
plt3.plot(x_expe_25s/D, y_expe_25s/D, "ro", markersize=marker_size)

# Legend
plt1.legend(loc='upper left')

# Figure saving
plt.savefig('Figures/' + figname + '.png', dpi=100, bbox_inches='tight')
