#
# Import section
#
import subprocess
import sys
import numpy as np
import fluidfoam

from analytic_coulomb2D import analytic_coulomb2D


def RMS(YM, UM, Yex, Uex):
    UexI = np.interp(YM, Yex, Uex, right=np.nan)
    RMSU = np.nanstd(UM[:] - UexI[:])
    return RMSU


#
#
#
zmin = 0.0
zmax = 0.06
#
# compute the analytical solution in dimensionless form
#
nx = 120
xex = np.linspace(0, 1.0, nx)
# dimensionless parameters
mus = 0.24
phi0 = 0.595
eta_e = 1.0 + 2.5 * phi0
# dimensional parameters
D = 0.06
rho_f = 950.0
rho_p = 1050.0
drho = rho_p - rho_f
etaf = 2.105e-5 * rho_f
g = 9.81

hp = 0.045825 / D
# pressure gradient
dpdx = -80e0 / (drho * g)

# Compute the analytical solution
alphaex = np.ones(nx) * phi0
toto = np.where(xex[:] > hp)
alphaex[toto] = 0.0

pex = np.zeros(nx)
for i in range(nx):
    if alphaex[nx - i - 1] > 0.0:
        pex[nx - i - 1] = pex[nx - i] + alphaex[nx - i] * (
            xex[nx - i] - xex[nx - i - 1]
        )

[uex, hc] = analytic_coulomb2D(nx, xex, dpdx, hp, mus, phi0, eta_e)

duxmax = 0.0
nuex = np.zeros(nx)
for i in range(nx - 1):
    duexdz = (uex[i] - uex[i - 1]) / (xex[i] - xex[i - 1])
    duxmax = max([duxmax, duexdz])
    nuex[i] = mus * pex[i] / (rho_p * (np.abs(duexdz) + 1e-6))
#
# dimensional form
#
U0 = drho * g * D**2 / etaf
uex = uex * U0
xex = xex * D
pex = pex * drho * g * D

print("max(uex)=" + str(np.max(uex)) + " m/s")

#########################################
# Reading SedFoam results
#########################################
sol = "./"
proc = subprocess.Popen(
    ["foamListTimes", "-case", sol, "-latestTime"], stdout=subprocess.PIPE
)
output = proc.stdout.read()
tread = output.decode().rstrip()
if float(tread) > 1900:
    tread = "1900"
tread = tread + "/"
tread = "1900/"

X, Y, Z = fluidfoam.readmesh(sol)
alpha = fluidfoam.readscalar(sol, tread, "alpha.a")
Ua = fluidfoam.readvector(sol, tread, "U.a")
Ub = fluidfoam.readvector(sol, tread, "U.b")
# pff = fluidfoam.readscalar(sol, tread, 'pff')
# p = fluidfoam.readscalar(sol, tread, 'p')

Ny = np.size(Y)
U = np.zeros(Ny)
U = alpha[:] * Ua[0, :] + (1 - alpha[:]) * Ub[0, :]

print("max(Ub)=" + str(np.amax(Ub)) + " m/s")
# figure(1)
# plot(Ub[0,:],Y)
# plot(uex,xex)
# show()
RMSU = RMS(Y, Ub[0, :], xex, uex)
print("RMS U=", RMSU)
assert RMSU <= 3e-3

# =============================================================================
