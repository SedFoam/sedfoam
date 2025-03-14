from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os


case = "3DOscillSheetFlow"
basepath = "../LES/"

sol = basepath + case + "/"

# plot creation and parameters
matplotlib.rcParams.update({"font.size": 20})
matplotlib.rc("text", usetex=True)
plt.figure(figsize=[10, 12])
gs = plt.GridSpec(5, 4)
gs.update(left=0.08, right=0.975, top=0.95, bottom=0.15, wspace=0.2, hspace=1.0)
font_legend = matplotlib.font_manager.FontProperties(
    family="Comic Sans MS", weight="bold", style="normal", size=14
)

period = 5
omega = 2 * np.pi / period
nu = 1e-6
delta = np.sqrt(2 * nu / omega)
intbed = 0.0915 * 0.5 / 0.606 - 0.0265

ymin = -0.0115 / delta
ymax = 0.02 / delta
#
# Freestream velocity
#
# plot 1
ax1 = plt.subplot(gs[0, :])
ax1.axis([0, 360, -1.4, 1.4])
ax1.set_ylabel(r"$u^f_\infty/U^f_m$")
ax1.set_xlabel(r"Phase (deg.)")
ax1.grid(linestyle=":")

#
# Concentration
#
# plot 2
ax2 = plt.subplot(gs[1:3, 0])
ax2.axis([0, 0.65, ymin, ymax])
ax2.set_title("A")
ax2.set_xlabel(r"$\langle \bar\phi \rangle$")
ax2.set_ylabel(r"$y/\delta$")
ax2.grid(linestyle=":")

# plot 3
ax3 = plt.subplot(gs[1:3, 1])
ax3.axis([0, 0.65, ymin, ymax])
ax3.set_title("B")
ax3.set_xlabel(r"$\langle \bar\phi \rangle$")
ax3.set_yticklabels([])
ax3.grid(linestyle=":")

# plot 4
ax4 = plt.subplot(gs[1:3, 2])
ax4.axis([0, 0.65, ymin, ymax])
ax4.set_title("C")
ax4.set_xlabel(r"$\langle \bar\phi \rangle$")
ax4.set_yticklabels([])
ax4.grid(linestyle=":")

# plot 5
ax5 = plt.subplot(gs[1:3, 3])
ax5.axis([0, 0.65, ymin, ymax])
ax5.set_title("D")
ax5.set_xlabel(r"$\langle \bar\phi \rangle$")
ax5.set_yticklabels([])
ax5.grid(linestyle=":")

#
# Velocity
#
# plot 6
ax6 = plt.subplot(gs[3:5, 0])
ax6.axis([0, 1.6, ymin, ymax])
ax6.set_xlabel(r"$\langle\tilde u^f\rangle_F/U^f_m$")
ax6.set_ylabel(r"$y/\delta$")
ax6.grid(linestyle=":")

# plot 7
ax7 = plt.subplot(gs[3:5, 1])
ax7.axis([0, 1.6, ymin, ymax])
ax7.set_xlabel(r"$\langle\tilde u^f\rangle_F/U^f_m$")
ax7.set_yticklabels([])
ax7.grid(linestyle=":")

# plot 8
ax8 = plt.subplot(gs[3:5, 2])
ax8.axis([0, 1.6, ymin, ymax])
ax8.set_xlabel(r"$\langle\tilde u^f\rangle_F/U^f_m$")
ax8.set_yticklabels([])
ax8.grid(linestyle=":")

# plot 9
ax9 = plt.subplot(gs[3:5, 3])
ax9.axis([0, 1.6, ymin, ymax])
ax9.set_xlabel(r"$\langle\tilde u^f\rangle_F/U^f_m$")
ax9.set_yticklabels([])
ax9.grid(linestyle=":")


# Plot freestream velocity
period = 5
omega = 2 * np.pi / period
nu = 1e-6
delta = np.sqrt(2 * nu / omega)
t = np.linspace(0, 5, 101)
sin = np.sin(omega * t)
ax1.plot(t * 360 / period, sin, "k")
ax1.plot(t[0] * 360 / period, sin[0], "ro", markersize=10, markeredgecolor="maroon")
ax1.plot(t[13] * 360 / period, sin[13], "ro", markersize=10, markeredgecolor="maroon")
ax1.plot(t[25] * 360 / period, sin[25], "ro", markersize=10, markeredgecolor="maroon")
ax1.plot(t[37] * 360 / period, sin[37], "ro", markersize=10, markeredgecolor="maroon")
ax1.text(6, -0.8, r"A")
ax1.text(44, -0.5, r"B")
ax1.text(88, -0.3, r"C")
ax1.text(132, -0.5, r"D")

# Numerical data reading
avg_file_phi = Dataset(os.path.join(sol, "postProcessing", case + "_alpha.nc"))
avg_file_ub = Dataset(os.path.join(sol, "postProcessing", case + "_Ub.nc"))

phi = avg_file_phi.variables["alpha"][:]
umean_f = avg_file_ub.variables["Ub"][:]
y = avg_file_phi.variables["pos"][:]

# Plot data
# Concentration
ax2.plot(phi[:, 0], (y - intbed) / delta, "-k")
ax3.plot(phi[:, 3], (y - intbed) / delta, "-k")
ax4.plot(phi[:, 5], (y - intbed) / delta, "-k")
ax5.plot(phi[:, 8], (y - intbed) / delta, "-k")

# Velocity
ax6.plot(umean_f[0, :, 0], (y - intbed) / delta, "-k")
ax7.plot(umean_f[0, :, 3], (y - intbed) / delta, "-k")
ax8.plot(umean_f[0, :, 5], (y - intbed) / delta, "-k")
ax9.plot(umean_f[0, :, 8], (y - intbed) / delta, "-k")

plt.savefig(
    "Figures/tuto3DOscillSheetFlow.png", facecolor="w", edgecolor="w", format="png"
)
plt.show()
