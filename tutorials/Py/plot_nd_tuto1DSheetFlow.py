import subprocess
import sys
import numpy as np
import fluidfoam
from pylab import figure, subplot, grid, axis, xlabel, ylabel, show
from pylab import xticks, savefig
from pylab import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib import rc
import matplotlib as mpl


def readOpenFoam(sol):
    import subprocess

    #
    # Reading SedFoam results
    try:
        proc = subprocess.Popen(
            ["foamListTimes", "-latestTime", "-case", sol], stdout=subprocess.PIPE
        )
    except FileNotFoundError:
        print("foamListTimes : command not found")
        print("Do you have load OpenFoam environement?")
        sys.exit(0)
    output = proc.stdout.read()
    tread = output.decode().rstrip().split("\n")[0]
    nt = 1
    X, Y, Z = fluidfoam.readmesh(sol)
    alpha = fluidfoam.readscalar(sol, tread, "alpha.a")
    Ua = fluidfoam.readvector(sol, tread, "U.a")
    Ub = fluidfoam.readvector(sol, tread, "U.b")
    Tauf = fluidfoam.readtensor(sol, tread, "Taub")
    Taus = fluidfoam.readtensor(sol, tread, "Taua")
    k = fluidfoam.readscalar(sol, tread, "k.b")
    Theta = fluidfoam.readscalar(sol, tread, "Theta")

    return nt, Y, Ua[0, :], Ub[0, :], alpha, Tauf[3, :], Taus[3, :], k, Theta


#
# Change fontsize
#
rc("text", usetex=True)
#
# linewidth
#
lw = 3
ms = 8
matplotlib.rcParams.update({"font.size": 20})
mpl.rcParams["lines.linewidth"] = lw
mpl.rcParams["lines.markersize"] = ms

#
gs = gridspec.GridSpec(1, 4)
gs.update(left=0.065, right=0.975, top=0.95, bottom=0.15, wspace=0.15, hspace=0.2)
#
# Figure size
#
figwidth = 14
figheight = 6

#
# Parameters
#
rho_f = 1000.0
rho_s = 1190.0
rs = 1.5e-3
alphasmax = 0.65

Ustarexp = 0.05
#
# load experimental results
#
[Zexpp, Phiexp] = np.loadtxt("DATA/Revil-Baudard2015_phi.txt", skiprows=1, unpack=True)
[Zexp, Uexp, uwexp] = np.loadtxt(
    "DATA/Revil-Baudard2015_uRxz.txt", skiprows=1, unpack=True
)

#########################################
#
# Loading OpenFoam results
#
basepath = "../RAS/"

Nx = 1
Ny = 400
Nz = 1

casedir = "1DSheetFlow/"
label = "K-Epsilon + Kinetic Theory"

############################################
# case 1
#

solpath = basepath + casedir
nt, y, ua, ub, alpha0, Tauxy, Tausxy, k, Theta = readOpenFoam(solpath)

#
# Plots
#

zmin = -2
zmax = 25

Z0num = -0.01 - rs

umax = 1.0/0.05
#
#
# figure 1
#
#
fig = figure(num=1, figsize=(figwidth, figheight), dpi=60, facecolor="w", edgecolor="w")

x_ticks = np.array([0, 1, 2, 3])
x_labels = ["0", "1", "2", "3"]
#
# ax1
#
ax1 = subplot(gs[0, 0])
p1_5 = ax1.plot(ub, (y*0.17 + Z0num) / (2.0 * rs)*0.003/0.17, "-b")
p1_4 = ax1.plot(ua, (y*0.17 + Z0num) / (2.0 * rs)*0.003/0.17, "-r")
p1_1 = ax1.plot(Uexp[:]/0.05, Zexp[:] / (2.0 * rs)*0.003/0.17, "ok", markersize=ms)
ylabel(r"$y/d_p$")
xlabel(r"$U/u_\tau$")
axis([0, umax, zmin*0.003/0.17, zmax*0.003/0.17])
grid()
#
# ax2
#
ax2 = subplot(gs[0, 1])
p2_3 = ax2.plot(alpha0, (y*0.17 + Z0num) / (2.0 * rs), "-r", label=label)
p2_1 = ax2.plot(Phiexp, Zexpp / (2.0 * rs), "ok", markersize=ms, label="experiment")
ax2.set_xscale("log")
handles2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(handles2[0:4], labels2[0:4], loc=1, prop={"size": 14})
ax2.set_yticklabels([])
xlabel(r"$\alpha$")
axis([1e-4, 1, zmin, zmax])
grid()

#
# ax3
#
ax3 = subplot(gs[0, 2])
p3_2 = ax3.plot(Tauxy, (y*0.17 + Z0num) / (2.0 * rs), "-b")
p3_2 = ax3.plot(Tausxy, (y*0.17 + Z0num) / (2.0 * rs), "--r")
p3_1 = ax3.plot(
    rho_f * uwexp/1000/0.05**2, Zexp / (2.0 * rs), "ok", markersize=ms, label="experiment"
)
ax3.set_yticklabels([])
xlabel(r"$\tau_{xz}^\ast$")
xticks(x_ticks, x_labels)
axis([0, 3/1000/0.05**2, zmin, zmax])
grid()

#
# ax4
#
ax4 = subplot(gs[0, 3])
p4_2 = ax4.plot(k, (y*0.17 + Z0num) / (2.0 * rs), "-b", label="TKE: k")
p4_2 = ax4.plot(
    1.5 * Theta, (y*0.17 + Z0num) / (2.0 * rs), "--r", label=r"Gran. Temp.: $3/2\theta$"
)
# p4_1 = ax4.plot(rho_f * uwexp, Zexp / (2. * rs), 'ok',
#                markersize=ms, label='experiment')
ax4.legend(fontsize=12)
ax4.set_yticklabels([])
xlabel(r"$k^\ast, 3/2\theta^\ast$")
axis([0, 7.5e-3/0.05**2, zmin, zmax])
grid()

savefig("Figures/res1_nd_tuto3.png", facecolor="w", edgecolor="w", format="png")

show()
