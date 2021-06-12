import subprocess
import sys
import numpy as np
import fluidfoam
import matplotlib.pyplot as plt
plt.ion()

############### Plot properties #####################
import matplotlib.ticker as mticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
 
label_size = 20
legend_size = 12
fontsize=25
linewidth=2
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['legend.fontsize'] = legend_size
plt.rcParams['lines.linewidth'] = linewidth
plt.rcParams['axes.labelsize'] = fontsize
####################################################

######################
# Load DEM data
######################

zDEM, phiDEM, vxPDEM, vxFDEM = np.loadtxt('DATA/BedloadTurbDEM.txt', unpack=True)

######################
#Read SedFoam results
######################
sol = '../1DBedLoadTurb/'
try:
    proc = subprocess.Popen(
            ["foamListTimes", "-latestTime", "-case", sol],
            stdout=subprocess.PIPE,)
except:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read() #to obtain the output of function foamListTimes from the subprocess
timeStep = output.decode().rstrip().split('\n')[0] #Some management on the output to obtain a number

#Read the data
X, Y, Z = fluidfoam.readmesh(sol)
z = Y
phi = fluidfoam.readscalar(sol, timeStep, 'alpha_a')
vxPart = fluidfoam.readvector(sol, timeStep, 'Ua')[0]
vxFluid = fluidfoam.readvector(sol, timeStep, 'Ub')[0]

######################
#Plot results
######################
d = 0.006 #6mm diameter particles
plt.figure(figsize=[10,5])
plt.subplot(141)
plt.plot(phiDEM, zDEM/d, 'k--', label=r'DEM')
plt.plot(phi, z/d, label=r'SedFoam')
plt.xlabel(r'$\phi$', fontsize=25)
plt.ylabel(r'$\frac{z}{d}$', fontsize=30, rotation=True, horizontalalignment='right')
plt.grid()
plt.legend()

plt.subplot(142)
I = np.where(phiDEM>0.001)[0]
plt.plot(vxPDEM[I], zDEM[I]/d, 'k--', label=r'DEM')
I = np.where(phi>0.001)[0]
plt.plot(vxPart[I], z[I]/d, label=r'SedFoam')
plt.xlabel(r'$v_x^p$', fontsize=25)
plt.ylim([-1.525, 32.025])
plt.grid()
ax = plt.gca()
ax.set_yticklabels([])

plt.subplot(143)
plt.plot(phiDEM*vxPDEM, zDEM/d, 'k--', label=r'DEM')
plt.plot(phi*vxPart, z/d, label=r'SedFoam')
plt.xlabel(r'$q = \phi v_x^p$', fontsize=25)
plt.grid()
ax = plt.gca()
ax.set_yticklabels([])

plt.subplot(144)
plt.plot(vxFDEM, zDEM/d, 'k--', label=r'DEM')
plt.plot(vxFluid, z/d, label=r'SedFoam')
plt.xlabel(r'$v_x^f$', fontsize=25)
plt.grid()
ax = plt.gca()
ax.set_yticklabels([])

plt.savefig('Figures/res_TutoBedloadTurb.png', bbox_inches='tight')
