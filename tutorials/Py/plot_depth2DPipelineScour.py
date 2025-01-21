import numpy as np
from pylab import mlab
import matplotlib.pyplot as plt
import fluidfoam
from scipy.interpolate import griddata


# Function reading the mesh dimensions
def readxy(case):
    x, y, z = fluidfoam.readmesh(case)
    return x, y


# Function returning the bed profile
def max_depth(case, t0, tfinal, dt, xi, yi):
    ybed = np.zeros(ngridx)
    arr_size = int((tfinal - t0) / dt + 0.01)
    yscour = np.zeros(arr_size)
    time = np.zeros(arr_size)
    i = -1
    j = 10
    t = t0
    while t < tfinal:
        t = t + dt
        i = i + 1
        if np.mod(t, 1) == 0:
            timename = str(int(t)) + "/"
        else:
            timename = str(t) + "/"
        a = fluidfoam.readscalar(case, timename, "alpha.a")
        ai = griddata((x, y), a, (xi, yi), method="linear")
        for j in range(ngridx):
            ybed[j] = np.max(yi[(np.where(ai[:, j] > 0.5))])
        yscour[i] = np.min(ybed[:])
        time[i] = t
    return yscour, time


# Case information
case = "../RAS/2DPipelineScour/"

# Number of division for linear interpolation
ngridx = 500
ngridy = 180

# Interpolation grid dimensions
xinterpmin = -0.05
xinterpmax = 0.2
yinterpmin = -0.075
yinterpmax = 0.015

# Interpolation grid
xi = np.linspace(xinterpmin, xinterpmax, ngridx)
yi = np.linspace(yinterpmin, yinterpmax, ngridy)
xinterp, yinterp = np.meshgrid(xi, yi)

# Maximum depth calculation
x, y = readxy(case)
t0 = 0
tfinal = 30
dt = 0.5
max_depth, time = max_depth(case, t0, tfinal, dt, xinterp, yinterp)

# Experimental data collection
expe = "DATA/Mao_depth_expe.txt"
time_expe, max_depth_expe = np.loadtxt(expe, usecols=(0, 1), unpack=True, delimiter=",")

# -------------PLOT SECTION------------- #
# Figure dimensions and parameters
fig_sizex = 18
fig_sizey = 9
font_size = 40
figname = "maximum_depth"
line_style = "-"
line_color = "C1"
line_width = 4
label_num = "sedFoam"
marker_size = 16
label_expe = "Mao [1986]"

# Figure creation
fig = plt.figure(figsize=(fig_sizex, fig_sizey), dpi=100)
plt.rcParams.update({"font.size": font_size})

# Plot axis
time_min = 0
time_max = 35
ymin = 0
ymax = 1
plt.axis([time_min, time_max, ymin, ymax])
plt.xlabel("t (s)")
plt.ylabel("S/D")
plt.grid()

# Numerical results plotting
plt.plot(
    time[0:60],
    abs((max_depth[0:60] + 0.025) / 0.05),
    label=label_num,
    linewidth=line_width,
    ls=line_style,
    color=line_color,
)

# Experimental data plotting
plt.plot(time_expe, max_depth_expe, "ro", label=label_expe, markersize=marker_size)

# Legend
plt.legend(loc="lower right", frameon=False, prop={"size": 28})

# Figure saving
plt.savefig("Figures/" + figname + ".png", dpi=100, bbox_inches="tight")
