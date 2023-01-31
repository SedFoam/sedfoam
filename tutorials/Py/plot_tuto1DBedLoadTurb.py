import subprocess
import sys
import numpy as np
import fluidfoam
import matplotlib.pyplot as plt

#             Plot properties
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from matplotlib import rc

plt.ion()

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc("text", usetex=True)

label_size = 20
legend_size = 12
fontsize = 25
linewidth = 2
plt.rcParams["xtick.labelsize"] = label_size
plt.rcParams["ytick.labelsize"] = label_size
plt.rcParams["legend.fontsize"] = legend_size
plt.rcParams["lines.linewidth"] = linewidth
plt.rcParams["axes.labelsize"] = fontsize
####################################################

######################
# Load DEM data
######################

zDEM, phiDEM, vxPDEM, vxFDEM, TDEM = np.loadtxt("DATA/BedloadTurbDEM.txt", unpack=True)

######################
# Read SedFoam results
######################
sol = "../RAS/1DBedLoadTurb/"
try:
    proc = subprocess.Popen(
        ["foamListTimes", "-latestTime", "-case", sol],
        stdout=subprocess.PIPE,
    )
except FileNotFoundError:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)

# to obtain the output of function foamListTimes from the subprocess
output = proc.stdout.read()

# Some management on the output to obtain a number
timeStep = output.decode().rstrip().split("\n")[0]

# Read the data
X, Y, Z = fluidfoam.readmesh(sol)
z = Y
phi = fluidfoam.readscalar(sol, timeStep, "alpha.a")
vxPart = fluidfoam.readvector(sol, timeStep, "U.a")[0]
vxFluid = fluidfoam.readvector(sol, timeStep, "U.b")[0]
T = fluidfoam.readscalar(sol, timeStep, "Theta")

######################
# Plot results
######################
d = 0.006  # 6mm diameter particles
plt.figure(figsize=[10, 5])
plt.subplot(141)
plt.plot(phiDEM, zDEM / d, "k--", label=r"DEM")
plt.plot(phi, z / d, label=r"SedFoam")
plt.xlabel(r"$\phi$", fontsize=25)
plt.ylabel(r"$\frac{z}{d}$", fontsize=30, rotation=True, horizontalalignment="right")
plt.grid()
plt.ylim([-1.525, 32.025])
plt.legend()

plt.subplot(142)
IndPhi = np.where(phiDEM > 0.001)[0]
plt.plot(vxPDEM[IndPhi], zDEM[IndPhi] / d, "r--")
IndPhi = np.where(phi > 0.001)[0]
plt.plot(vxPart[IndPhi], z[IndPhi] / d, "r", label=r"$v_x^p$")
plt.plot(vxFDEM, zDEM / d, "b--")
plt.plot(vxFluid, z / d, "b", label=r"$u_x^f$")
plt.xlabel(r"$v_x^p$, $u_x^f$", fontsize=25)
plt.ylim([-1.525, 32.025])
plt.grid()
plt.legend()
ax = plt.gca()
ax.set_yticklabels([])
plt.legend()

plt.subplot(143)
plt.plot(phiDEM * vxPDEM, zDEM / d, "k--", label=r"DEM")
plt.plot(phi * vxPart, z / d, label=r"SedFoam")
plt.xlabel(r"$q = \phi v_x^p$", fontsize=25)
plt.grid()
plt.ylim([-1.525, 32.025])
ax = plt.gca()
ax.set_yticklabels([])

plt.subplot(144)

IndPhi = np.where(phiDEM > 0.001)[0]
plt.plot(TDEM[IndPhi], zDEM[IndPhi] / d, "k--", label=r"DEM")
IndPhi = np.where(phi > 0.001)[0]
plt.plot(T[IndPhi], z[IndPhi] / d, label=r"SedFoam")
plt.xlabel(r"$T$", fontsize=25)
plt.grid()
plt.ylim([-1.525, 32.025])
ax = plt.gca()
ax.set_yticklabels([])

plt.savefig("Figures/res_TutoBedloadTurb.png", bbox_inches="tight")

plt.show(block=True)
