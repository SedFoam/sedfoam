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
except:
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

phi_interp = np.interp(zDATA, z, phi)
rms_phi = rms(phi_interp - phiDATA)
assert rms_phi <= 0.02

vxPart_interp = np.interp(zDATA, z, vxPart)
I = np.where(zDATA < 17.5 * 0.006)
rms_vxP = rms(vxPart_interp[I] - vxPDATA[I])
assert rms_vxP <= 0.1

vxFluid_interp = np.interp(zDATA, z, vxFluid)
rms_vxF = rms(vxFluid_interp - vxFDATA)
assert rms_vxF <= 0.15

T_interp = np.interp(zDATA, z, T)
I = np.where(zDATA < 17.5 * 0.006)
rms_T = rms(T_interp[I] - TDATA[I])
assert rms_T <= 0.02
