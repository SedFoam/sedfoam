import subprocess
import sys
import numpy as np
import fluidfoam
from pylab import plt, matplotlib, show


plt.rc("text", usetex=True)
font = {"family": "serif", "size": 16, "serif": ["computer modern roman"]}
plt.rc("font", **font)
plt.rc("legend", **{"fontsize": 16})

##################################################################
# Parameters to retreive dimensionless variables
##################################################################

d = 160e-6  # particle diameter in m
gravity = 9.81  # gravity in m/s2
rhoFluid = 1041  # fluid density in kg/m3
h = 0.0049  # initial granular height in m

timeAdim = (d / gravity) ** 0.5
velAdim = 1000.0 * (gravity * d) ** 0.5
pressureAdim = rhoFluid * h * gravity


##################################################################
# Experimental data extracted from Pailha et al. (2008)
##################################################################

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/v_t_exp_562.txt", delimiter="\t", names=True
)
time_v_562 = data1["time"] / timeAdim
velocity_562 = data1["velocity"] / velAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/p_t_exp_562.txt", delimiter="\t", names=True
)
time_p_562 = data1["time"] / timeAdim
pressure_562 = data1["pressure"] / pressureAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/v_t_exp_568.txt", delimiter="\t", names=True
)
time_v_568 = data1["time"] / timeAdim
velocity_568 = data1["velocity"] / velAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/p_t_exp_568.txt", delimiter="\t", names=True
)
time_p_568 = data1["time"] / timeAdim
pressure_568 = data1["pressure"] / pressureAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/v_t_exp_578.txt", delimiter="\t", names=True
)
time_v_578 = data1["time"] / timeAdim
velocity_578 = data1["velocity"] / velAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/p_t_exp_578.txt", delimiter="\t", names=True
)
time_p_578 = data1["time"] / timeAdim
pressure_578 = data1["pressure"] / pressureAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/v_t_exp_584.txt", delimiter="\t", names=True
)
time_v_584 = data1["time"] / timeAdim
velocity_584 = data1["velocity"] / velAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/p_t_exp_584.txt", delimiter="\t", names=True
)
time_p_584 = data1["time"] / timeAdim
pressure_584 = data1["pressure"] / pressureAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/v_t_exp_592.txt", delimiter="\t", names=True
)
time_v_592 = data1["time"] / timeAdim
velocity_592 = data1["velocity"] / velAdim

data1 = np.genfromtxt(
    "DATA/ExperimentalDataPailha2008/p_t_exp_592.txt", delimiter="\t", names=True
)
time_p_592 = data1["time"] / timeAdim
pressure_592 = data1["pressure"] / pressureAdim


#########################################
# Loading SedFoam results
#########################################
sol = "../laminar/1DWetAvalanche"

try:
    proc = subprocess.Popen(
        ["foamListTimes", "-latestTime", "-case", sol], stdout=subprocess.PIPE
    )
except FileNotFoundError:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
final_tread = output.decode().rstrip().split("\n")[0]


X, Y, Z = fluidfoam.readmesh(sol)
time_sim_dila_0 = []
vel_sim_dila_0 = []
p_sim_dila_0 = []
tolAlpha = 0.54

for i in range(200, int(final_tread)):
    if i % 20 == 0:
        time_sim_dila_0.append((i - 200) / timeAdim)
        tread = str(i) + "/"
        alpha_0 = fluidfoam.readscalar(sol, tread, "alpha.a")
        Ua_0 = fluidfoam.readvector(sol, tread, "U.a")
        p_rbgh_0 = fluidfoam.readscalar(sol, tread, "p_rbgh")
        for k in range(len(alpha_0)):
            if alpha_0[k] < tolAlpha:
                vel_sim_dila_0.append(Ua_0[0, k] * 1000 / velAdim)
                p_sim_dila_0.append(p_rbgh_0[4] / pressureAdim)
                break


#########################################
# 				Plots
#########################################


# time - velocity plot including the experimental data (Pailha et al 2008) and the numerical results

plt.figure()
plt.plot(
    time_v_562,
    velocity_562,
    marker="o",
    markersize=0,
    linestyle="--",
    color="lightpink",
    label="Exp. - $\\phi_0 = 0.562$",
)
plt.plot(
    time_v_568,
    velocity_568,
    marker="o",
    markersize=0,
    linestyle="--",
    color="plum",
    label="Exp. - $\\phi_0 = 0.568$",
)
plt.plot(
    time_v_578,
    velocity_578,
    marker="o",
    markersize=0,
    linestyle="--",
    color="mediumpurple",
    label="Exp. - $\\phi_0 = 0.578$",
)
plt.plot(
    time_v_584,
    velocity_584,
    marker="o",
    markersize=0,
    linestyle="-",
    color="royalblue",
    label="Exp. - $\\phi_0 = 0.584$",
)
plt.plot(
    time_v_592,
    velocity_592,
    marker="o",
    markersize=0,
    linestyle="-",
    color="navy",
    label="Exp. - $\\phi_0 = 0.592$",
)

plt.plot(
    time_sim_dila_0,
    vel_sim_dila_0,
    marker="o",
    markersize=5,
    linestyle="",
    linewidth=1.2,
    color="r",
    label="SedFoam  - $\\phi_0 = 0.592$",
)

plt.ylabel("$\\frac{v^s}{\\sqrt{gd}}$ [$-$]", fontsize=18)
plt.xlabel("$\\frac{t}{\\sqrt{d/g}}$ [$-$]", fontsize=18)
plt.axis([-1000, 160000, -0.00005, 0.03001])
plt.grid()
plt.tight_layout()
plt.savefig("Figures/velocityPlot1D_phi0592" + ".png", dpi=200)


# time - pressure plot including the experimental data (Pailha et al 2008) and the numerical results

plt.figure()
plt.plot(
    time_p_562,
    pressure_562,
    marker="o",
    markersize=0,
    linestyle="--",
    color="lightpink",
    label="Exp. - $\\phi_0 = 0.562$",
)
plt.plot(
    time_p_568,
    pressure_568,
    marker="o",
    markersize=0,
    linestyle="--",
    color="plum",
    label="Exp. - $\\phi_0 = 0.568$",
)
plt.plot(
    time_p_578,
    pressure_578,
    marker="o",
    markersize=0,
    linestyle="--",
    color="mediumpurple",
    label="Exp. - $\\phi_0 = 0.578$",
)
plt.plot(
    time_p_584,
    pressure_584,
    marker="o",
    markersize=0,
    linestyle="-",
    color="royalblue",
    label="Exp. - $\\phi_0 = 0.584$",
)
plt.plot(
    time_p_592,
    pressure_592,
    marker="o",
    markersize=0,
    linestyle="-",
    color="navy",
    label="Exp. - $\\phi_0 = 0.592$",
)

plt.plot(
    time_sim_dila_0,
    p_sim_dila_0,
    marker="o",
    markersize=5,
    linestyle="",
    linewidth=1.2,
    color="r",
    label="SedFoam  - $\\phi_0 = 0.592$",
)

plt.ylabel("$\\frac{p^f}{\\rho^f g h_o}$ [$-$]", fontsize=18)
plt.xlabel("$\\frac{t}{\\sqrt{d/g}}$ [$-$]", fontsize=18)
plt.axis([-1000, 160000, -0.12, 0.12])
plt.grid()
plt.tight_layout()
plt.legend(prop={"size": 12.0}, loc=0)
plt.savefig("Figures/pressurePlot1D_phi0592.png", dpi=200)


show(block=True)
