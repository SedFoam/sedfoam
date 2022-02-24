import subprocess
import sys
import numpy as np
import fluidfoam
from pylab import matplotlib, mpl, figure, subplot, savefig, show
import matplotlib.gridspec as gridspec
from analytic_coulomb2D import analytic_coulomb2D

#
# Change fontsize
#
matplotlib.rcParams.update({'font.size': 20})
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['lines.markeredgewidth'] = 1
#
# Change subplot sizes
#
gs = gridspec.GridSpec(1, 3)
gs.update(left=0.1, right=0.95, top=0.95,
          bottom=0.1, wspace=0.125, hspace=0.25)
#
# Figure size
#
figwidth = 18
figheight = 9
#
#
#
zmin = 0.
zmax = 0.065
#
# compute the analytical solution in dimensionless form
#
nx = 200
xex = np.linspace(0, 1., nx)
# dimensionless parameters
mus = 0.32
phi0 = 0.6
eta_e = (1. + 2.5 * phi0)
# dimensional parameters
D = 0.065
etaf = 0.27
rho_f = 1070.
rho_p = 1190.
drho = rho_p - rho_f
g = 9.81

hp = 0.03225 / D
print("hp=" + str(hp))
# pressure gradient
dpdx = -100. / (drho * g)

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
#########################################
# Loading OpenFoam results
#########################################
#
case = '1DBedLoad'
basepath = '../laminar/'
# basepath='../../'
sol = basepath + case + '/'

try:
    proc = subprocess.Popen(
        ['foamListTimes', '-latestTime', '-case', sol], stdout=subprocess.PIPE)
except:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
tread = output.decode().rstrip().split('\n')[0]

Nx = 1
Ny = 200
Nz = 1

eps_file = sol + case + '.eps'

#########################################
# Reading SedFoam results
#########################################

X, Y, Z = fluidfoam.readmesh(sol)
alpha = fluidfoam.readscalar(sol, tread, 'alpha.a')
Ua = fluidfoam.readvector(sol, tread, 'U.a')
Ub = fluidfoam.readvector(sol, tread, 'U.b')
pff = fluidfoam.readscalar(sol, tread, 'pff')
tau = fluidfoam.readtensor(sol, tread, 'Taua')[3]
#p = fluidfoam.readscalar(sol, tread, 'p')

Ny = np.size(Y)
U = np.zeros(Ny)
U = alpha[:] * Ua[0, :] + (1 - alpha[:]) * Ub[0, :]

print("max(Ub)=" + str(np.amax(Ub)) + " m/s")

#########################################
# figure 1
#########################################

figure(num=1, figsize=(figwidth, figheight),
       dpi=60, facecolor='w', edgecolor='w')

ax1 = subplot(gs[0, 0])
l11, = ax1.plot(alpha[:], Y[:], '-r')
l1, = ax1.plot(alphaex[:], xex[:], '--k')
ax1.set_ylabel('y (m)')
ax1.set_xlabel(r'$\alpha$')
ax1.set_xlim(0,  np.max(np.max(alpha)) * 1.1)
ax1.set_ylim(zmin, zmax)

ax2 = subplot(gs[0, 1])
l21, = ax2.plot(U[:], Y[:], '-r')
l2, = ax2.plot(uex[:], xex[:], '--k')
ax2.set_xlabel('u ($m/s$)')
ax2.set_xlim(0,  np.max(uex) * 1.1)
ax2.set_ylim(zmin, zmax)
ax2.set_yticklabels([''])

ax3 = subplot(gs[0, 2])
l31, = ax3.plot(pff[:], Y[:], '-r')
l3, = ax3.plot(pex[:], xex[:], '--k')
ax3.set_xlabel('p ($N/m^2$)')
ax3.set_xlim(0,  np.max(pex) * 1.1)
ax3.set_ylim(zmin, zmax)
ax3.set_yticklabels([''])

savefig('Figures/res1_tuto2.png', facecolor='w', edgecolor='w', format='png')

show(block=True)

# Fix Python 2.x.
try: input = raw_input
except NameError: pass
toto = input("Hit a key to close the figure")
