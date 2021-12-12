import subprocess
import numpy as np
import fluidfoam
from pylab import *
import matplotlib.gridspec as gridspec

def rms(x):
    return np.sqrt(x.dot(x)/x.size)

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
          bottom=0.2, wspace=0.125, hspace=0.25)
#
# Figure size
#
figwidth = 15
figheight = 6
#
#
#
zmin = 0
zmax = 1
#
# compute the analytical solution in dimensionless form
#
# dimensional parameters
d = 3.175
rho_f = 1e-3
rho_p = 1.
drho = rho_p - rho_f
g = 0.0

[phiexp,yexpphi] =np.genfromtxt('DATA/simpleShearVescovi2014_phiY.dat',delimiter=';',comments='#', unpack=True)
[uexp,yexp] =np.genfromtxt('DATA/simpleShearVescovi2014_uY.dat',delimiter=';',comments='#', unpack=True)
[thetaexp,yexptheta] =np.genfromtxt('DATA/simpleShearVescovi2014_TY.dat',delimiter=';',comments='#', unpack=True)



######################################
# Loading OpenFoam results
#########################################
case = '1DSimpleShear'
basepath = '../laminar/'
sol = basepath + case + '/'

#
# Reading SedFoam results
#
try:
    proc = subprocess.Popen(
        ['foamListTimes', '-latestTime', '-case', sol], stdout=subprocess.PIPE)
except:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
tread = output.decode().rstrip().split('\n')[0]

eps_file = sol + case + '.eps'

#########################################
# Reading SedFoam results
#########################################
prec=9
X,Y,Z = fluidfoam.readmesh(sol,True,precision=prec)
nx,ny,nz = X.shape

alpha = fluidfoam.readscalar(sol, tread, 'alpha_a',True,precision=prec)
Ua    = fluidfoam.readvector(sol, tread, 'Ua',True,precision=prec)
pa    = fluidfoam.readscalar(sol, tread, 'pa',True,precision=prec)
Theta = fluidfoam.readscalar(sol, tread, 'Theta',True,precision=prec)

Ny = np.size(Y)
H = np.max(np.max(Y))
U = 1

#########################################
# figure 1
#########################################

figure(num=1, figsize=(figwidth, figheight),
       dpi=60, facecolor='w', edgecolor='w')

iprof=0
ax1 = subplot(gs[0, 0])
l11, = ax1.plot(alpha[0,:,iprof], Y[0,:,iprof]/H, '-r')
l12, = ax1.plot(phiexp, yexpphi, 'ob',label='simple shear')
l13, = ax1.plot([0.57, 0.57], [0,1], '--k')

ax1.set_ylabel('y/H [-]')
ax1.set_xlabel(r'$\alpha$')
ax1.set_xlim(0,  np.max(np.max(alpha)) * 1.1)
ax1.set_ylim(zmin, zmax)

ax2 = subplot(gs[0, 1])
l21, = ax2.plot(Ua[0,0,:,iprof]/U, Y[0,:,iprof]/H, '-r',label='sedFoam')
l22, = ax2.plot(uexp, yexp, 'ob',label='simple shear')
ax2.legend(loc='upper left',fontsize=12)
ax2.set_xlabel(r'$u/U$')
ax2.set_xlim(0,  1.)
ax2.set_ylim(zmin, zmax)
ax2.set_yticklabels([''])

ax3 = subplot(gs[0, 2])
l31, = ax3.semilogx(Theta[0,:,iprof]/U, Y[0,:,iprof]/H, '-r')
l32, = ax3.plot(thetaexp, yexptheta, 'ob')
ax3.set_xlabel(r'$\theta/U$')
ax3.set_xlim(1e-4,  1e-2)
ax3.set_ylim(zmin, zmax)
ax3.set_yticklabels([''])


#savefig('Figures/res1_DryAvalanche.png', facecolor='w', edgecolor='w', format='png')
show(block=True)

phi_interp = np.interp(Y[0,:,iprof]/H,yexpphi, phiexp);
rms_phi = rms(phi_interp - alpha[0,:,iprof])
assert(rms_phi<=0.03)

u_interp = np.interp(Y[0,:,iprof]/H,yexp, uexp);
rms_u = rms(u_interp - Ua[0,0,:,iprof])
assert(rms_u<=0.04)

theta_interp = np.interp(Y[0,:,iprof]/H,yexptheta, thetaexp);
rms_theta = rms(theta_interp - Theta[0,:,iprof])/np.mean(thetaexp)
assert(rms_theta<=0.4)



