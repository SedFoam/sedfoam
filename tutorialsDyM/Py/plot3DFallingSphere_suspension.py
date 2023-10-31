import numpy as np
from pylab import matplotlib, plt, show
import os
import sys

DynamycMesh="Morphing"
#DynamycMesh="Overset"
if DynamycMesh=="Overset":
	DynamycMeshV=DynamycMesh+"/sphereAndBackground"
else:
        DynamycMeshV=DynamycMesh
D=0.015 #diameter of the sphere

#################### loading experimental data ##############
data_TenCate_displacement = np.genfromtxt("DATA/tenCate2002displacementRe1.5.csv", delimiter=",", names=["X", "Y"])
data_TenCate_vel = np.genfromtxt("DATA/tenCate2002velocityRe1.5.csv", delimiter=",", names=["X", "Y"])

######## loading overSedDymFoam data without sediments using overset ########
time_noSediment = np.genfromtxt("DATA/fallingSphereNoSedimentOverset/time.txt",delimiter='\t', names=True)
v_noSediment = np.genfromtxt("DATA/fallingSphereNoSedimentOverset/vz.txt",delimiter='\t', names=True)
z_noSediment = np.genfromtxt("DATA/fallingSphereNoSedimentOverset/zcenter.txt",delimiter='\t', names=True)

time_noSediment_v=[]
v_noSediment_v=[]
z_noSediment_v=[]

for i in range(10,len(time_noSediment)-10):
	time_noSediment_v.append(time_noSediment[i][0])
	v_noSediment_v.append(v_noSediment[i][0])
	z_noSediment_v.append(z_noSediment[i][0]/D)
	
######## loading overSedDymFoam data without sediments using morphing mesh########
time_noSediment_morph_d = np.genfromtxt("DATA/fallingSphereNoSedimentMorphing/time.txt",delimiter='\t', names=True)
v_noSediment_morph_d = np.genfromtxt("DATA/fallingSphereNoSedimentMorphing/vz.txt",delimiter='\t', names=True)
z_noSediment_morph_d = np.genfromtxt("DATA/fallingSphereNoSedimentMorphing/zcenter.txt",delimiter='\t', names=True)

time_noSediment_morph=[]
v_noSediment_morph=[]
z_noSediment_morph=[]

for i in range(10,len(time_noSediment_morph_d)-10):
	time_noSediment_morph.append(time_noSediment_morph_d[i][0])
	v_noSediment_morph.append(v_noSediment_morph_d[i][0])
	z_noSediment_morph.append(z_noSediment_morph_d[i][0]/D)

############## loading overSedDymFoam results ###########
sol="../FallingSphereSuspension"+DynamycMeshV+'/'
if os.path.exists(sol+'time')==False:
        print ("To generate the postprocessing files you need to execute makeFiles in the case first")   
        sys.exit()
datatime = np.genfromtxt(sol+'time',delimiter='\t', names=True)
dataVz= np.genfromtxt(sol+'vz',delimiter='\t', names=True)
dataPosz = np.genfromtxt(sol+'zcenter',delimiter='\t', names=True)


time_Sim=[]
v_Sim=[]
z_Simt=[]

count=0
for i in range(10,len(dataVz)-10):
	if i%1==0:
		v_Sim.append(dataVz[i][0])
		time_Sim.append(datatime[i][0])
		z_Simt.append(dataPosz[i][0]/D)

######################### plots ################################

plt.rcParams.update({'font.size': 30})
plt.rc('font', family='serif')
plt.rc('text',usetex=True)
font = {'family':'serif','size':30, 'serif': ['computer modern roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':30})
matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']
matplotlib.rcParams["legend.framealpha"] = None

limX=0.7

fig, axs = plt.subplots(2, 1, figsize=(16.4, 16.4))
#axs[0].set_title('Sphere position')
axs[0].plot(time_noSediment_v, z_noSediment_v,color='b', linestyle='--')
axs[0].plot(time_noSediment_morph, z_noSediment_morph,color='cyan', linestyle='--')
axs[0].plot(time_Sim, z_Simt,color='brown', linestyle='--')
axs[0].plot(data_TenCate_displacement['X'], data_TenCate_displacement['Y']-8, linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[0].set_xlim([-0, limX])
axs[0].set_ylim([ -1.51,0])
axs[0].set_xticklabels([])
axs[0].set_ylabel('z/D [-]',fontsize="30")
axs[0].grid()

#axs[1].set_title('Sphere velocity')
axs[1].plot(time_noSediment_v, v_noSediment_v,color='b', linestyle='--',label="overSedDyMFoam - OversetMesh - $\\phi=0.0$")
axs[1].plot(time_noSediment_morph, v_noSediment_morph,color='cyan', linestyle='--',label="overSedDyMFoam - MorphingMesh - $\\phi=0.0$")
axs[1].plot(time_Sim, v_Sim,color='brown', linestyle='--',label="overSedDyMFoam - %s - $\\phi=0.05$"%DynamycMesh)
axs[1].plot(data_TenCate_vel['X'], data_TenCate_vel['Y'], linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[1].set_xlabel('t [s]',fontsize="30")
axs[1].set_ylabel('Vertical velocity [m/s]',fontsize="30")
axs[1].set_xlim([-0, limX])
axs[1].set_ylim([ -0.0401,0.0001])
axs[1].grid()
axs[1].legend(loc="upper right",fontsize="25")
fig.tight_layout()
fig.savefig('Figures/FallingSphereSuspension.png')
plt.show()



