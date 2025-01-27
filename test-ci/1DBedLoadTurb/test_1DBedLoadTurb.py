import subprocess
import sys
import numpy as np
import fluidfoam


def rms(x):
    return np.sqrt(x.dot(x) / x.size)


##################
# Loading Data
##################
zDATA, phiDATA, vxPDATA, vxFDATA, TDATA = np.loadtxt(
    "DATA/BedloadTurbDATA.txt", unpack=True
)

######################################
# Loading OpenFoam results
#########################################
case = "./"
basepath = "./"
sol = basepath + case + "/"

#
# Reading SedFoam results
#
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

#########################################
# Reading SedFoam results
#########################################
X, Y, Z = fluidfoam.readmesh(sol)

z = Y
phi = fluidfoam.readscalar(sol, tread, "alpha.a")
vxPart = fluidfoam.readvector(sol, tread, "U.a")[0]
vxFluid = fluidfoam.readvector(sol, tread, "U.b")[0]
T = fluidfoam.readscalar(sol, tread, "Theta")

rms_phi = rms(phi - phiDATA)
assert rms_phi <= 0.02

Ind = np.where(phi > 1e-4)
rms_vxP = rms(vxPart[Ind] - vxPDATA[Ind])
assert rms_vxP <= 0.1

rms_vxF = rms(vxFluid - vxFDATA)
assert rms_vxF <= 0.1

rms_T = rms(T[Ind] - TDATA[Ind])
assert rms_T <= 0.02
