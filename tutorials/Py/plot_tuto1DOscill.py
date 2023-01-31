#
# Import section
#
import numpy as np
from pylab import figure, subplot, show, matplotlib, ylabel, title, xlabel
from math import sqrt, pi
import fluidfoam

matplotlib.rcParams.update({"font.size": 12})

#########################################
#
# Physical parameters
#
h = 0.0451
U0 = 1.1
T0 = 4.0
viscof = 1e-6

delta = sqrt(viscof * T0 / pi)

#
# Loading OpenFoam results
#
mypath = "../RAS/"

#
#  Solution folder
#
basepath = "1DOscill/"

############################################
# case 1
#
sol = mypath + basepath

x, y, z = fluidfoam.readmesh(sol)

datalist = ["1.5", "2", "2.5", "3"]
datalist = ["3.5", "4", "4.5", "5"]
datalist = ["7.5", "8", "8.5", "9"]
# datalist=['11.5' , '12','12.5','13']
# datalist=['15.5' , '16','16.5','17']

#
# Loading literature results
#

#
# U
#
npzfiles = np.load("DATA/G03_U_fig3a.npz")
#
# DNS model
#
u_DNS_phi0 = npzfiles["arr_0"]
zu_DNS_phi0 = npzfiles["arr_1"]
u_DNS_phi315 = npzfiles["arr_2"]
zu_DNS_phi315 = npzfiles["arr_3"]
u_DNS_phi45 = npzfiles["arr_4"]
zu_DNS_phi45 = npzfiles["arr_5"]
u_DNS_phi90 = npzfiles["arr_6"]
zu_DNS_phi90 = npzfiles["arr_7"]

#
# Original Wilcox k-w model
#
u_OW_phi0 = npzfiles["arr_8"]
zu_OW_phi0 = npzfiles["arr_9"]
u_OW_phi315 = npzfiles["arr_10"]
zu_OW_phi315 = npzfiles["arr_11"]
u_OW_phi45 = npzfiles["arr_12"]
zu_OW_phi45 = npzfiles["arr_13"]
u_OW_phi90 = npzfiles["arr_14"]
zu_OW_phi90 = npzfiles["arr_15"]
#
# Guizien k-w model
#
u_phi0 = npzfiles["arr_16"]
zu_phi0 = npzfiles["arr_17"]
u_phi315 = npzfiles["arr_18"]
zu_phi315 = npzfiles["arr_19"]
u_phi45 = npzfiles["arr_20"]
zu_phi45 = npzfiles["arr_21"]
u_phi90 = npzfiles["arr_22"]
zu_phi90 = npzfiles["arr_23"]

#
# k
#
npzfiles = np.load("DATA/G03_k_fig5.npz")
#
# DNS model
#
TKE_DNS_phi0 = npzfiles["arr_0"]
zTKE_DNS_phi0 = npzfiles["arr_1"]
TKE_DNS_phi315 = npzfiles["arr_2"]
zTKE_DNS_phi315 = npzfiles["arr_3"]
TKE_DNS_phi45 = npzfiles["arr_4"]
zTKE_DNS_phi45 = npzfiles["arr_5"]
TKE_DNS_phi90 = npzfiles["arr_6"]
zTKE_DNS_phi90 = npzfiles["arr_7"]
#
# Original Wilcox k-w model
#
TKE_OW_phi0 = npzfiles["arr_8"]
zTKE_OW_phi0 = npzfiles["arr_9"]
TKE_OW_phi315 = npzfiles["arr_10"]
zTKE_OW_phi315 = npzfiles["arr_11"]
TKE_OW_phi45 = npzfiles["arr_12"]
zTKE_OW_phi45 = npzfiles["arr_13"]
TKE_OW_phi90 = npzfiles["arr_14"]
zTKE_OW_phi90 = npzfiles["arr_15"]
#
# Guizien k-w model
#
TKE_phi0 = npzfiles["arr_16"]
zTKE_phi0 = npzfiles["arr_17"]
TKE_phi315 = npzfiles["arr_18"]
zTKE_phi315 = npzfiles["arr_19"]
TKE_phi45 = npzfiles["arr_20"]
zTKE_phi45 = npzfiles["arr_21"]
TKE_phi90 = npzfiles["arr_22"]
zTKE_phi90 = npzfiles["arr_23"]

#
# uw
#
npzfiles = np.load("DATA/G03_uw_fig4.npz")
#
# DNS model
#
R_DNS_phi0 = npzfiles["arr_0"]
zR_DNS_phi0 = npzfiles["arr_1"]
R_DNS_phi315 = npzfiles["arr_2"]
zR_DNS_phi315 = npzfiles["arr_3"]
R_DNS_phi45 = npzfiles["arr_4"]
zR_DNS_phi45 = npzfiles["arr_5"]
R_DNS_phi90 = npzfiles["arr_6"]
zR_DNS_phi90 = npzfiles["arr_7"]
#
# Original Wilcox k-w model
#
R_OW_phi0 = npzfiles["arr_8"]
zR_OW_phi0 = npzfiles["arr_9"]
R_OW_phi315 = npzfiles["arr_10"]
zR_OW_phi315 = npzfiles["arr_11"]
R_OW_phi45 = npzfiles["arr_12"]
zR_OW_phi45 = npzfiles["arr_13"]
R_OW_phi90 = npzfiles["arr_14"]
zR_OW_phi90 = npzfiles["arr_15"]
#
# Guizien k-w model
#
R_phi0 = npzfiles["arr_16"]
zR_phi0 = npzfiles["arr_17"]
R_phi315 = npzfiles["arr_18"]
zR_phi315 = npzfiles["arr_19"]
R_phi45 = npzfiles["arr_20"]
zR_phi45 = npzfiles["arr_21"]
R_phi90 = npzfiles["arr_22"]
zR_phi90 = npzfiles["arr_23"]

phi = ["-45", "0", "45", "90"]

#%%  graphes a tracer

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
    zu = y - np.min(y)

    figure(1)

    ax1 = subplot(1, 1, 1)
    pO = ax1.plot(u / U0, zu / delta, "--r", label="OpenFOAM")  # phi="+phi[i])

    if i == 0:
        p1 = ax1.plot(
            u_phi315, zu_phi315, "-b", label="Guizien et al (2003)"
        )  # phi=-45")
        p2 = ax1.plot(u_OW_phi315, zu_OW_phi315, "-.g", label="Wilcox")  # phi=-45")
        p3 = ax1.plot(u_DNS_phi315, zu_DNS_phi315, ":k", label="DNS")  # phi=-45")

        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles, labels, loc=1)
    elif i == 1:
        p1 = ax1.plot(u_phi0, zu_phi0, "-b", label="Guizien et al (2003) phi=0")
        p2 = ax1.plot(u_OW_phi0, zu_OW_phi0, "-.g", label="Wilcox phi=0")
        p3 = ax1.plot(u_DNS_phi0, zu_DNS_phi0, ":k", label="DNS phi=0")
    elif i == 2:
        p1 = ax1.plot(u_phi45, zu_phi45, "-b", label="Guizien et al (2003) phi=45")
        p2 = ax1.plot(u_OW_phi45, zu_OW_phi45, "-.g", label="Wilcox phi=45")
        p3 = ax1.plot(u_DNS_phi45, zu_DNS_phi45, ":k", label="DNS phi=45")
    elif i == 3:
        p1 = ax1.plot(u_phi90, zu_phi90, "-b", label="Guizien et al (2003) phi=90")
        p2 = ax1.plot(u_OW_phi90, zu_OW_phi90, "-.g", label="Wilcox phi=90")
        p3 = ax1.plot(u_DNS_phi90, zu_DNS_phi90, ":k", label="DNS phi=90")

    ylabel(r"$z / \delta$")
    ax1.axis([-1.2, 1.2, 0, 30])

#
#
# TKE,R
#
#
i = -1
sub = [1, 3, 2, 4]
for data in datalist:
    i = i + 1
    print(data)

    k = fluidfoam.readscalar(sol, data + "/", "k.b")
    zk = y - np.min(y)
    # omega = fluidfoam.readscalar(sol, data+'/', 'omega')

    #
    #
    #
    figure(2)

    ax2 = subplot(2, 2, sub[i])
    pO = ax2.plot(k / U0**2, zk / delta, "--r", label="OpenFOAM")
    # pO   = ax2.plot(omega/U0**3,zk/delta,'-r',label="OpenFOAM")

    if i == 0:
        p1 = ax2.plot(TKE_phi90, zTKE_phi90, "-b", label="Guizien et al (2003)")
        p2 = ax2.plot(TKE_OW_phi90, zTKE_OW_phi90, "-.g", label="Wilcox")
        p3 = ax2.plot(TKE_DNS_phi90, zTKE_DNS_phi90, ":k", label="DNS")
        ax2.axis([0, 0.008, 0, 30])
        title("phi=-45")
        ylabel(r"$z / \delta$")
        ax2.set_xticklabels([""])

    elif i == 1:
        p1 = ax2.plot(TKE_phi315, zTKE_phi315, "-b", label="Guizien et al (2003)")
        p2 = ax2.plot(TKE_OW_phi315, zTKE_OW_phi315, "-.g", label="Wilcox")
        p3 = ax2.plot(TKE_DNS_phi315, zTKE_DNS_phi315, ":k", label="DNS")
        ax2.axis([0, 0.003, 0, 30])
        title("phi=0")
        xlabel(r"$k / U_0^2$")
        ylabel(r"$z / \delta$")
    elif i == 2:
        p1 = ax2.plot(TKE_phi0, zTKE_phi0, "-b", label="Guizien et al (2003)")
        p2 = ax2.plot(TKE_OW_phi0, zTKE_OW_phi0, "-.g", label="Wilcox")
        p3 = ax2.plot(TKE_DNS_phi0, zTKE_DNS_phi0, ":k", label="DNS")
        ax2.axis([0, 0.009, 0, 30])
        title("phi=45")
        ax2.set_yticklabels([""])
        ax2.set_xticklabels([""])

    elif i == 3:
        p1 = ax2.plot(TKE_phi45, zTKE_phi45, "-b", label="Guizien et al (2003)")
        p2 = ax2.plot(TKE_OW_phi45, zTKE_OW_phi45, "-.g", label="Wilcox")
        p3 = ax2.plot(TKE_DNS_phi45, zTKE_DNS_phi45, ":k", label="DNS")
        ax2.axis([0, 0.016, 0, 30])
        title("phi=90")
        ax2.set_yticklabels([""])
        xlabel(r"$k / U_0^2$")

        # handles, labels = ax2.get_legend_handles_labels()
        # ax2.legend(handles, labels,loc=1)

    ax2.axis([0, 0.015, 0, 40])

    #
    # Ruw
    #
    figure(3)

    Taub = fluidfoam.readtensor(sol, data + "/", "Taub")
    Rfw = -Taub[3, :] / 1e3
    zR = y - np.min(y)

    ax3 = subplot(2, 2, sub[i])
    pO = ax3.plot(Rfw / U0**2, zR / delta, "--r", label="OpenFOAM")

    if i == 0:
        p1 = ax3.plot(R_phi315, zR_phi315, "-b", label="Guizien et al (2003)")
        p2 = ax3.plot(R_OW_phi315, zR_OW_phi315, "-.g", label="Wilcox")
        p3 = ax3.plot(R_DNS_phi315, zR_DNS_phi315, ":k", label="DNS")
        title("phi=-45")
        ylabel(r"$z / \delta$")
        ax3.set_xticklabels([""])

    elif i == 1:
        p1 = ax3.plot(R_phi0, zR_phi0, "-b", label="Guizien et al (2003)")
        p2 = ax3.plot(R_OW_phi0, zR_OW_phi0, "-.g", label="Wilcox")
        p3 = ax3.plot(R_DNS_phi0, zR_DNS_phi0, ":k", label="DNS")
        title("phi=0")
        xlabel(r"$Ruw / U_0^2$")
        ylabel(r"$z / \delta$")
    elif i == 2:
        p1 = ax3.plot(R_phi45, zR_phi45, "-b", label="Guizien et al (2003)")
        p2 = ax3.plot(R_OW_phi45, zR_OW_phi45, "-.g", label="Wilcox")
        p3 = ax3.plot(R_DNS_phi45, zR_DNS_phi45, ":k", label="DNS")
        title("phi=45")
        ax3.set_yticklabels([""])
        ax3.set_xticklabels([""])

    elif i == 3:
        p1 = ax3.plot(R_phi90, zR_phi90, "-b", label="Guizien et al (2003)")
        p2 = ax3.plot(R_OW_phi90, zR_OW_phi90, "-.g", label="Wilcox")
        p3 = ax3.plot(R_DNS_phi90, zR_DNS_phi90, ":k", label="DNS")
        title("phi=90")
        ax3.set_yticklabels([""])
        xlabel(r"$Ruw / U_0^2$")

    # handles, labels = ax3.get_legend_handles_labels()
    # ax3.legend(handles, labels,loc=1)

    ax3.axis([-2.5e-3, 2.5e-3, 0, 40])


show()
