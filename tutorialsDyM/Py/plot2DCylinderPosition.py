import numpy as np
from fluidfoam import readmesh, readvector, readscalar, readtensor
from pylab import matplotlib, plt, show
import csv

#Particles
diameterPart = 1.5e-3 #Diameter of the particles composing the bed, in m
diameterPartI = 3.3*diameterPart   #Diameter of the intruder, in m

#################### loading DEM results  ##############

runA = np.genfromtxt("DATA/Yade_A.csv", delimiter=",", names=["time", "posZ"])
runB = np.genfromtxt("DATA/Yade_B.csv", delimiter=",", names=["time", "posZ"])
runC = np.genfromtxt("DATA/Yade_C.csv", delimiter=",", names=["time", "posZ"])
runD = np.genfromtxt("DATA/Yade_D.csv", delimiter=",", names=["time", "posZ"])

#################### loading overSedDymFoam results  ##############

sol = '../2DCylinderUniformGranularFlow/'

datatime= np.genfromtxt(sol+'time',delimiter='\t', names=True)
dataPosz = np.genfromtxt(sol+'zcenter',delimiter='\t', names=True)

timeSed=[]
Z_sed=[]
for i in range(10,len(dataPosz)-10):
	timeSed.append(datatime[i][0])
	Z_sed.append( dataPosz[i][0]/diameterPartI)
		
##################### plot ####################

plt.rcParams.update({'font.size': 16})
plt.rc('font', family='serif')
plt.rc('text',usetex=True)
font = {'family':'serif','size':16, 'serif': ['computer modern roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':16})
matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']
matplotlib.rcParams["legend.framealpha"] = None

plt.rcParams["figure.figsize"] = (4,6)
plt.xlabel('Time [s]',fontsize="16")
plt.ylabel('z/D [-]',fontsize="16")
plt.plot(runA["time"],runA["posZ"],color='orange',alpha=0.25,linestyle='-')
plt.plot(runB["time"],runB["posZ"],color='orange',alpha=0.5,linestyle='-' )
plt.plot(runC["time"],runC["posZ"],color='orange',alpha=0.75,linestyle='-' ,label='YADE' )
plt.plot(runD["time"],runD["posZ"],color='orange',alpha=0.99,linestyle='-' )
plt.plot(timeSed, Z_sed,color='b',alpha=0.99,linestyle='--',label='overSedDymFoam' )

plt.grid()
plt.xlim([-0.05,8])
plt.ylim([-0.5,2])
plt.legend(loc="upper left",fontsize=12)
plt.tight_layout()
plt.savefig('Figures/CylinderGranularFlow_PositionVSTime.png', dpi=600)
plt.show()
