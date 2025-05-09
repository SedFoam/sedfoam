import numpy as np
from pylab import figure, subplot, xlabel, ylabel, show
from math import sqrt, pi
import fluidfoam

#########################################
#
# Physical parameters
#
h = 0.01
U0 = 1.1
T0 = 4.0
viscof = 1e-6

delta = sqrt(viscof * T0 / pi)

# Loading OpenFoam results
#
mypath = "../laminar/"
#
#  Solution folder
#
basepath = "1DStokesBL/"

############################################
# case 1
#
sol = mypath + basepath

x, y, z = fluidfoam.readmesh(sol)

datalist = ["3.5", "4", "4.5", "5", "5.5"]

#       graphs

#
#
# U
#
#
i = -1

zex = np.linspace(0, h, 200)


for data in datalist:
    i = i + 1
    print(data)

    # ufile=basepath+data+'Uf.xy'
    # zu,u=np.loadtxt(ufile,unpack=True)
    Ub = fluidfoam.readvector(sol, data + "/", "U.b")
    u = Ub[0, :]
    Taub = fluidfoam.readtensor(sol, data + "/", "Taub")
    Rfw = Taub[3, :] / 1e3
    z = y - np.min(y)

    tex = float(data)
    uex = U0 * (
        np.sin(2 * np.pi * tex / T0)
        - np.exp(-zex / delta) * np.sin(2 * np.pi * tex / T0 - zex / delta)
    )
    tauex = (
        viscof
        * U0
        / delta
        * np.exp(-zex / delta)
        * (
            np.sin(2 * np.pi * tex / T0 - zex / delta)
            + np.cos(2 * np.pi * tex / T0 - zex / delta)
        )
    )

    figure(1)

    ax1 = subplot(1, 1, 1)
    pO = ax1.plot(u / U0, z / delta, "--r", label="OpenFOAM")  # phi="+phi[i])
    p = ax1.plot(uex / U0, zex / delta, "-m", label="Stokes")  # phi="+phi[i])

    if i == 0:
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles, labels, loc=1)

    xlabel(r"$u^f / U_0$")
    ylabel(r"$z / \delta$")

    ax1.axis([-1.2, 1.2, 0, 10])

    # Ruw
    #
    figure(2)

    ax3 = subplot(1, 1, 1)
    pO = ax3.plot(Rfw / U0**2, z / delta, "--r", label="OpenFOAM")
    p = ax3.plot(tauex / U0**2, zex / delta, "-m", label="Stokes")  # phi="+phi[i])
    xlabel(r"$\tau^f / \rho U_0^2$")
    ylabel(r"$z / \delta$")
    ax3.axis([-1.25e-3, 1.25e-3, 0, 10])


show()
