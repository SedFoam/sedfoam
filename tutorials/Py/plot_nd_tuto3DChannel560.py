from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os


case = "3DChannel560"
basepath = "../LES/"
sol = basepath + case + "/"

# plot creation and parameters
matplotlib.rcParams.update({"font.size": 20})
matplotlib.rc("text", usetex=True)
plt.figure(figsize=[15, 6])
gs = plt.GridSpec(1, 3)
gs.update(left=0.08, right=0.975, top=0.95, bottom=0.15, wspace=0.1, hspace=0.1)
font_legend = matplotlib.font_manager.FontProperties(
    family="Comic Sans MS", weight="bold", style="normal", size=14
)

h = 0.02
retau=560
nu=1e-6
utau=retau*nu/h

# plot 1
ax1 = plt.subplot(gs[0, 0])
ax1.axis([1e-10, 0.3, 0, 1])
ax1.set_xlabel(r"$\langle\bar\phi\rangle$")
ax1.set_ylabel(r"$y/h$")
ax1.grid(linestyle=":")


# plot 2
ax2 = plt.subplot(gs[0, 1])
#ax2.axis([0, 0.65, 0, 1])
ax2.set_xlabel(r"$\langle\tilde u^f\rangle_F$, $\langle\tilde u^s\rangle_F$")
ax2.set_yticklabels([])
ax2.grid(linestyle=":")


# plot 3
ax3 = plt.subplot(gs[0, 2])
#ax3.axis([0, 1e-2, 0, 1])
ax3.axis([-0.35, 1.0, 0, 1])
ax3.set_xlabel(
    r"$\langle\tilde u^{f'} \tilde v^{f'}\rangle_F$, "
    r"$\langle\tilde u^{s'} \tilde v^{s'}\rangle_F$"
)
ax3.set_yticklabels([])
ax3.grid(linestyle=":")


# Numerical data reading
avg_file = Dataset(os.path.join(sol, "postProcessing", case + "_averaged.nc"))
phi = avg_file.variables["alpha.a"][:]
ndumean_s = avg_file.variables["uaf.a"][:]
ndumean_f = avg_file.variables["ubf.a"][:]
nduprim_s = avg_file.variables["uaprimf.a"][:]
nduprim_f = avg_file.variables["ubprimf.a"][:]
ndy = avg_file.variables["y"][:]


# Load the data from the file
data = np.loadtxt('./DATA/kasagi_re=650_cc_dns_umean.dat', skiprows=1)
x_dns_umean = data[:, 1]
y_dns_umean = data[:, 3]
data = np.loadtxt('./DATA/kasagi_re=650_cc_dns_restress.dat', skiprows=1)
x_dns_restress = data[:, 1]
y_dns_restress = data[:, 3]


# Plot data
ax1.semilogx(phi, ndy , "-k")
ax2.plot(ndumean_s[0], ndy , "-r", label="LES-Re=560-solid")
ax2.plot(ndumean_f[0], ndy , "-k", label="LES-Re=560-fluid")
ax2.plot(y_dns_umean, x_dns_umean, linestyle='--', color='b', label="DNS-Re=650-single-phase")
ax3.plot(-nduprim_s[3, :], ndy, "-r", label="LES-Re=560-solid")
ax3.plot(-nduprim_f[3, :], ndy , "-k", label="LES-Re=560-fluid")
ax3.plot(y_dns_restress, x_dns_restress, linestyle='--', color='b', label="DNS-Re=650-single-phase")


ax2.set_xlim(0, 22)
ax3.set_xlim(0, 1) 
ax2.legend(fontsize='18')
ax3.legend(fontsize='18')

plt.savefig("Figures/tuto3DChannel560_nd.png", facecolor="w", edgecolor="w", format="png")
plt.show()


