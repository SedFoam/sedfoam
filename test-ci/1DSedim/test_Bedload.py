#
# Import section
#
import subprocess
import sys
import numpy as np
import fluidfoam
from pylab import figure, subplot, axis, xlabel, ylabel, show, savefig, plot
from pylab import title, matplotlib
import matplotlib.gridspec as gridspec
import matplotlib as mpl
sys.path.append("../../tutorials/Py/") 

from analytic_coulomb2D import *

def RMS(YM,UM,Yex,Uex):

    UexI=np.interp(YM,Yex,Uex,right=np.NAN);

    RMSU = np.nanstd(UM[:]-UexI[:])

    return RMSU

#
#
#
zmin = 0.
zmax = 0.06
#
# compute the analytical solution in dimensionless form
#
nx = 120
xex = np.linspace(0, 1., nx)
# dimensionless parameters
mus = 0.24
phi0 = 0.595
eta_e = (1. + 2.5 * phi0)
# dimensional parameters
D = 0.06
rho_f = 950.
rho_p = 1050.
drho = rho_p - rho_f
etaf = 2.105e-5*rho_f
g = 9.81

hp = 0.045825/D 
# pressure gradient
dpdx = -80e0 / (drho * g)

# Compute the analytical solution
alphaex = np.ones(nx) * phi0
toto = np.where(xex[:] > hp)
alphaex[toto] = 0.

pex = np.zeros(nx)
for i in range(nx):
    if alphaex[nx - i - 1] > 0.:
        pex[nx - i - 1] = pex[nx - i] + alphaex[nx - i] * \
            (xex[nx - i] - xex[nx - i - 1])

[uex, hc] = analytic_coulomb2D(nx, xex, dpdx, hp, mus, phi0, eta_e)

duxmax = 0.
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

#
# Change subplot sizes
#
gs = gridspec.GridSpec(1, 3)
gs.update(left=0.1, right=0.95, top=0.95,bottom=0.1, wspace=0.125, hspace=0.25)
#########################################
# Reading SedFoam results
#########################################
sol='./'
proc = subprocess.Popen(
    ['foamListTimes', '-case', sol, '-latestTime'], stdout=subprocess.PIPE)
output = proc.stdout.read()
tread = output.decode().rstrip()
if float(tread)>1900:
    tread='1900'
tread=tread+'/'
tread = '1900/'

X, Y, Z = fluidfoam.readmesh(sol)
alpha = fluidfoam.readscalar(sol, tread, 'alpha_a')
Ua = fluidfoam.readvector(sol, tread, 'Ua')
Ub = fluidfoam.readvector(sol, tread, 'Ub')
#pff = fluidfoam.readscalar(sol, tread, 'pff')
#p = fluidfoam.readscalar(sol, tread, 'p')

Ny = np.size(Y)
U = np.zeros(Ny)
U = alpha[:] * Ua[0, :] + (1 - alpha[:]) * Ub[0, :]

print("max(Ub)=" + str(np.amax(Ub)) + " m/s")
#figure(1)
#plot(Ub[0,:],Y)
#plot(uex,xex)
#show()
RMSU = RMS(Y,Ub[0,:],xex,uex)
print('RMS U=',RMSU)
assert(RMSU<=3e-3)

# =============================================================================
# #########################################
# # figure 1
# #########################################
# #
# # Change subplot sizes
# #
# gs = gridspec.GridSpec(1, 3)
# gs.update(left=0.1, right=0.95, top=0.95,
#           bottom=0.1, wspace=0.125, hspace=0.25)
# #
# # Figure size
# #
# figwidth = 18
# figheight = 9
# figure(num=1, figsize=(figwidth, figheight),
#        dpi=60, facecolor='w', edgecolor='w')
# 
# ax1 = subplot(gs[0, 0])
# l11, = ax1.plot(alpha[:], Y[:], '-r')
# l1, = ax1.plot(alphaex[:], xex[:], '--k')
# ax1.set_ylabel('y (m)')
# ax1.set_xlabel(r'$\alpha$')
# ax1.set_xlim(0,  np.max(np.max(alpha)) * 1.1)
# ax1.set_ylim(zmin, zmax)
# 
# ax2 = subplot(gs[0, 1])
# l21, = ax2.plot(U[:], Y[:], '-r')
# l2, = ax2.plot(uex[:], xex[:], '--k')
# ax2.set_xlabel('u ($m/s$)')
# ax2.set_xlim(0,  np.max(uex) * 1.1)
# ax2.set_ylim(zmin, zmax)
# ax2.set_yticklabels([''])
# 
# ax3 = subplot(gs[0, 2])
# l31, = ax3.plot(pff[:], Y[:], '-r')
# l3, = ax3.plot(pex[:], xex[:], '--k')
# ax3.set_xlabel('p ($N/m^2$)')
# ax3.set_xlim(0,  np.max(pex) * 1.1)
# ax3.set_ylim(zmin, zmax)
# ax3.set_yticklabels([''])
# 
# 
# show(block=True)
# =============================================================================

