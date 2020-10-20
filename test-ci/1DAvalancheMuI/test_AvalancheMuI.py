import analyticBagnold
import subprocess
import os
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
gs = gridspec.GridSpec(1, 7)
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
#
# compute the analytical solution in dimensionless form
#
# dimensional parameters
d = 0.02
rho_f = 1e-3
rho_p = 1
drho = rho_p - rho_f
g = 1

# dimensionless parameters
mus = 0.38
mu2 = 0.64
I0 = 0.3
Bphi = 0.31
phi0 = 0.59999999

beta = 26*np.pi/180.


#########################################
# Loading OpenFoam results
#########################################
case = './'
basepath = './'
sol = basepath + case + '/'

#
# Reading SedFoam results
#
try:
    proc = subprocess.Popen(
        ['foamListTimes', '-latestTime', '-case', sol], stdout=subprocess.PIPE)
except:
    print("foamListTimes : command not found")
    print("Did you loaded OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
tread = output.decode().rstrip().split('\n')[0]

#########################################
# Reading SedFoam results
#########################################

X,Y,Z = fluidfoam.readmesh(sol)
alpha = fluidfoam.readscalar(sol, tread, 'alpha_a')
Ua = fluidfoam.readvector(sol, tread, 'Ua')
Ub = fluidfoam.readvector(sol, tread, 'Ub')
pff = fluidfoam.readscalar(sol, tread, 'pff')
pa = fluidfoam.readscalar(sol, tread, 'pa')
p = fluidfoam.readscalar(sol, tread, 'p')
muI = fluidfoam.readscalar(sol, tread, 'muI')
nuEffa = fluidfoam.readscalar(sol, tread, 'nuEffa')
nuEffb = fluidfoam.readscalar(sol, tread, 'nuEffb')
nuFra = fluidfoam.readscalar(sol, tread, 'nuFra')
Tauf = fluidfoam.readtensor(sol, tread, 'Taub')
Taus = fluidfoam.readtensor(sol, tread, 'Taua')
try:
    gradUa = fluidfoam.readtensor(sol, tread, 'grad(Ua)')
except:
    print("grad(Ua) was not found -> Introduce - postProcess -func 'grad(Ua)' - in the command line")
    os.system("postProcess -func \'grad(Ua)\' -time "+tread)
    gradUa = fluidfoam.readtensor(sol, tread, 'grad(Ua)')


Ny = np.size(Y)
U = np.zeros(Ny)
U = alpha[:] * Ua[0, :] + (1 - alpha[:]) * Ub[0, :]

taua = np.zeros(Ny)
taub = np.zeros(Ny)
taua = Taus[3, :]
taub =  Tauf[3, :]
if (np.size(pff))==1:
    pff=pff*np.ones(Ny)

print("max(Ua)=" + str(np.amax(Ua)) + " m/s")

dudy = gradUa[3,:]
nu = np.zeros(Ny)
for i in range(Ny - 1):
    nu[i] = taua[0]/(dudy[i]*rho_p)

granularBed = np.where(pa+pff<(rho_p-rho_f)*g*d/10)
nx = np.size(granularBed)
H=Y[np.min(granularBed)]/d

print("H/d=",H)
zmax = np.max(Y)/d

[xex,I,alphaex,uex,pex,muIex,tauex,duexdz,nuex] = \
             analyticBagnold.analyticBagnold(nx,H,g,d,rho_p,rho_f,phi0,I0,Bphi,mus,mu2,beta)
print("max(uex)=" + str(np.max(uex)) + " m/s" + "max(phi)=",np.max(alphaex),"max(nuex)=",np.max(nuex))

#########################################
# figure 1
#########################################

# =============================================================================
# figure(num=1, figsize=(figwidth, figheight),
#        dpi=60, facecolor='w', edgecolor='w')
#
# ax1 = subplot(gs[0, 0])
# l11, = ax1.plot(alpha[:], Y[:]/d, '-r')
# l1, = ax1.plot(alphaex[:], xex[:]/d, '--k')
# ax1.set_ylabel('y/d [-]')
# ax1.set_xlabel(r'$\alpha$')
# ax1.set_xlim(0.5,  0.61)
# ax1.set_ylim(zmin, zmax)
#
# ax2 = subplot(gs[0, 1])
# l21, = ax2.plot(Ua[0,:]/(g*d)**0.5, Y[:]/d, '-r')
# l22, = ax2.plot(Ub[0,:]/(g*d)**0.5, Y[:]/d, '--b')
# l22, = ax2.plot((alpha*Ua[0,:]+(1-alpha)*Ub[0,:])/(g*d)**0.5, Y[:]/d, ':m')
# l2, = ax2.plot(uex[:]/(g*d)**0.5, xex[:]/d, '--k')
# 
# ax2.set_xlabel(r'$U/\sqrt{g d}$')
# ax2.set_xlim(0,  np.max([np.max(uex),np.max(Ua[0,:])])/(g*d)**0.5 * 1.1)
# ax2.set_ylim(zmin, zmax)
# ax2.set_yticklabels([''])
# 
# ax3 = subplot(gs[0, 2])
# l31, = ax3.plot((pff[:]+pa[:])/(rho_p*g*d), Y[:]/d, '-r')
# l31, = ax3.plot((pff[:])/(rho_p*g*d), Y[:]/d, ':r')
# l3, = ax3.plot(pex[:]/(rho_p*g*d), xex[:]/d, '--k')
# ax3.set_xlabel(r'$p/\rho_p g d$')
# ax3.set_xlim(0,  np.max(pex)/(rho_p*g*d) * 1.1)
# ax3.set_ylim(zmin, zmax)
# ax3.set_yticklabels([''])
# 
# ax4 = subplot(gs[0, 3])
# ax4.plot(muI,Y/d,'-r')
# ax4.plot(muIex,xex/d,'--k')
# ax4.plot([mus,mus],[0,H],':k')
# ax4.plot([mu2,mu2],[0,H],':k')
# ax4.set_xlabel(r'$\mu$(I) ')
# ax4.set_xlim(0.35,  0.7)
# ax4.set_ylim(zmin, zmax)
# ax4.set_yticklabels([''])
# 
# ax5 = subplot(gs[0, 4])
# l51, = ax5.plot(taua[:]/(rho_p*g*d), Y[:]/d, '-r')
# l52, = ax5.plot(taub[:]/(rho_p*g*d), Y[:]/d, '--b')
# l5, = ax5.plot(tauex[:]/(rho_p*g*d), xex[:]/d, '--k')
# ax5.set_xlabel(r'$\tau_s/\rho_p g d$')
# ax5.set_xlim(0,   max([np.max(taua),np.max(tauex)])/(rho_p*g*d) * 1.1)
# ax5.set_ylim(zmin, zmax)
# ax5.set_yticklabels([''])
# 
# ax6 = subplot(gs[0, 5])
# l61, = ax6.plot(nuFra[:]/(g*d**3)**0.5, Y[:]/d, '-r')
# l61, = ax6.plot(alpha[:]*nuEffa[0]/(g*d**3)**0.5, Y[:]/d, ':r')
# l61, = ax6.plot((1-alpha[:])*nuEffb[:]/(g*d**3)**0.5, Y[:]/d, ':b')
# l6, = ax6.plot(nuex[:]/(g*d**3)**0.5, xex[:]/d, '--k')
# ax6.set_xlabel(r'$\nu/\sqrt{g d^3}$')
# ax6.set_xlim(0,  max([np.max(nuex),np.max(nuFra)])/(g*d**3)**0.5 * 1.1)
# ax6.set_ylim(zmin, zmax)
# ax6.set_yticklabels([''])
# 
# ax7 = subplot(gs[0, 6])
# l71, = ax7.plot(dudy[:]/(g/d)**0.5, Y[:]/d, '-r',label='SedFoam')
# l7, = ax7.plot(duexdz[:]/(g/d)**0.5, xex[:]/d, '--k',label='Theoretical')
# ax7.set_xlabel(r'$\frac{d u^s}{dy}\sqrt{\frac{d}{g}}$ ')
# ax7.set_xlim(0,  np.max(duexdz)/(g/d)**0.5 * 1.1)
# ax7.set_ylim(zmin, zmax)
# ax7.set_yticklabels([''])
# ax7.legend(prop={'size':10.0},loc=0)

#show(block=True)

# =============================================================================

alpha_interp = np.interp(xex, Y[:],alpha);
rms_alpha = rms(alpha_interp - alphaex)
assert(rms_alpha<=0.05)

u_interp = np.interp(xex, Y[:],Ua[0,:]);
rms_u = rms(u_interp - uex)
assert(rms_u<=0.09)

muI_interp = np.interp(xex,Y[:], muI);
rms_muI = rms(muI_interp - muIex)
assert(rms_muI<=0.02)



