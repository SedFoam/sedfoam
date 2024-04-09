import numpy as np
from pylab import matplotlib, plt, show
import os
import sys

D=0.015 #diameter of the sphere

#################### loading experimental data ##############

data_TenCate_displacement = np.genfromtxt("DATA/tenCate2002displacementRe1.5.csv", delimiter=",", names=["X", "Y"])
data_TenCate_vel = np.genfromtxt("DATA/tenCate2002velocityRe1.5.csv", delimiter=",", names=["X", "Y"])

############## loading overSedDymFoam results ###########

sol="../FallingSphereMorphing/"
if os.path.exists(sol+'time.txt')==False:
        print ("b To generate the postprocessing files you need to execute makeFiles in the case first")   
        sys.exit()
datatime = np.genfromtxt(sol+'time.txt',delimiter='\t', names=True)
dataVz= np.genfromtxt(sol+'vz.txt',delimiter='\t', names=True)
dataPosz = np.genfromtxt(sol+'zcenter.txt',delimiter='\t', names=True)

time_Sim_morph=[]
v_Sim_morph=[]
z_Simt_morph=[]
count=0
for i in range(10,len(dataVz)-10):
	if i%1==0:
		v_Sim_morph.append(dataVz[i][0])
		time_Sim_morph.append(datatime[i][0])
		z_Simt_morph.append(dataPosz[i][0]/D)
        
sol="../FallingSphereSuspensionMorphing/"
if os.path.exists(sol+'time.txt')==False:
        print ("c To generate the postprocessing files you need to execute makeFiles in the case first")   
        sys.exit()
datatime = np.genfromtxt(sol+'time.txt',delimiter='\t', names=True)
dataVz= np.genfromtxt(sol+'vz.txt',delimiter='\t', names=True)
dataPosz = np.genfromtxt(sol+'zcenter.txt',delimiter='\t', names=True)

time_Sim_morph_sus=[]
v_Sim_morph_sus=[]
z_Simt_morph_sus=[]
count=0
for i in range(10,len(dataVz)-10):
	if i%1==0:
		v_Sim_morph_sus.append(dataVz[i][0])
		time_Sim_morph_sus.append(datatime[i][0])
		z_Simt_morph_sus.append(dataPosz[i][0]/D)

sol="../FallingSphereOverset/sphereAndBackground/"
if os.path.exists(sol+'time.txt')==False:
        print ("d To generate the postprocessing files you need to execute makeFiles in the case first")   
        sys.exit()
datatime = np.genfromtxt(sol+'time.txt',delimiter='\t', names=True)
dataVz= np.genfromtxt(sol+'vz.txt',delimiter='\t', names=True)
dataPosz = np.genfromtxt(sol+'zcenter.txt',delimiter='\t', names=True)

time_Sim_overset=[]
v_Sim_overset=[]
z_Simt_overset=[]
count=0
for i in range(10,len(dataVz)-10):
	if i%1==0:
		v_Sim_overset.append(dataVz[i][0])
		time_Sim_overset.append(datatime[i][0])
		z_Simt_overset.append(dataPosz[i][0]/D)

sol="../FallingSphereSuspensionOverset/sphereAndBackground/"
if os.path.exists(sol+'time.txt')==False:
        print ("e To generate the postprocessing files you need to execute makeFiles in the case first")   
        sys.exit()
datatime = np.genfromtxt(sol+'time.txt',delimiter='\t', names=True)
dataVz= np.genfromtxt(sol+'vz.txt',delimiter='\t', names=True)
dataPosz = np.genfromtxt(sol+'zcenter.txt',delimiter='\t', names=True)

time_Sim_overset_sus=[]
v_Sim_overset_sus=[]
z_Simt_overset_sus=[]
count=0
for i in range(10,len(dataVz)-10):
	if i%1==0:
		v_Sim_overset_sus.append(dataVz[i][0])
		time_Sim_overset_sus.append(datatime[i][0])
		z_Simt_overset_sus.append(dataPosz[i][0]/D)
        
######################### plots ################################

plt.rcParams.update({'font.size': 30})
plt.rc('font', family='serif')
plt.rc('text',usetex=True)
font = {'family':'serif','size':30, 'serif': ['computer modern roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':30})
matplotlib.rcParams["legend.framealpha"] = None

limX=0.7
alhpaV=0.6

fig, axs = plt.subplots(2, 2, figsize=(20, 16.4))
axs[0][0].text(0.25,0.05, 'Morphing mesh',  fontsize=30,rotation=0)
axs[0][0].plot(time_Sim_morph, z_Simt_morph,color='b',alpha=alhpaV, linestyle='-')
axs[0][0].plot(time_Sim_morph_sus, z_Simt_morph_sus,color='brown', linestyle='--')
axs[0][0].plot(data_TenCate_displacement['X'], data_TenCate_displacement['Y']-8, linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[0][0].set_xlim([-0, limX])
axs[0][0].set_ylim([ -1.51,0])
axs[0][0].set_xticklabels([])
axs[0][0].set_ylabel('z/D [-]',fontsize="30")
axs[0][0].grid()

axs[1][0].plot(time_Sim_morph, v_Sim_morph,color='b',alpha=alhpaV, linestyle='-',label="overSedDyMFoam - Morphing - $\\alpha=0.00$")
axs[1][0].plot(time_Sim_morph_sus, v_Sim_morph_sus,color='brown', linestyle='--',label="overSedDyMFoam - Morphing - $\\alpha=0.05$")
axs[1][0].plot(data_TenCate_vel['X'], data_TenCate_vel['Y'], linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[1][0].set_xlabel('t [s]',fontsize="30")
axs[1][0].set_ylabel('Vertical velocity [m/s]',fontsize="30")
axs[1][0].set_xlim([-0, limX])
axs[1][0].set_ylim([ -0.0401,0.0001])
axs[1][0].grid()
axs[1][0].legend(loc="upper right",fontsize="25")

axs[0][1].text(0.25,0.05, 'Overset mesh',  fontsize=30,rotation=0)
axs[0][1].plot(time_Sim_overset, z_Simt_overset,color='cyan', linestyle='-')
axs[0][1].plot(time_Sim_overset_sus, z_Simt_overset_sus,color='purple', linestyle='--')
axs[0][1].plot(data_TenCate_displacement['X'], data_TenCate_displacement['Y']-8, linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[0][1].set_xlim([-0, limX])
axs[0][1].set_ylim([ -1.51,0])
axs[0][1].set_xticklabels([])
axs[0][1].set_yticklabels([])
axs[0][1].grid()

axs[1][1].plot(time_Sim_overset, v_Sim_overset,color='cyan', linestyle='-',label="overSedDyMFoam - Overset - $\\alpha=0.00$")
axs[1][1].plot(time_Sim_overset_sus, v_Sim_overset_sus,color='purple', linestyle='--',label="overSedDyMFoam - Overset - $\\alpha=0.05$")
axs[1][1].plot(data_TenCate_vel['X'], data_TenCate_vel['Y'], linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[1][1].set_xlabel('t [s]',fontsize="30")
axs[1][1].set_xlim([-0, limX])
axs[1][1].set_ylim([ -0.0401,0.0001])
axs[1][1].set_yticklabels([])
axs[1][1].grid()
axs[1][1].legend(loc="upper right",fontsize="25")

fig.tight_layout()
fig.savefig('Figures/FallingSphereMorphAndOverset.png')
plt.show()
