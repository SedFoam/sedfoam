"""
Read and Plot a time series of OpenFoam postProcessing force
============================================================

This example reads and plots a series of postProcessing force
"""

###############################################################################
# Read the postProcessing files
# -----------------------------
#
# .. note:: In this example it reads and merges two postProcessing files
#           automatically (with the 'mergeTime' option)

# import readforce function from fluidfoam package
from fluidfoam.readpostpro import readforce
import numpy as np
#import probeLib as plib
import fluidfoam


H=0.015

data_TenCate_displacement = np.genfromtxt("DATA/tenCate2002displacementRe1.5.csv", delimiter=",", names=["X", "Y"])
data_TenCate_vel = np.genfromtxt("DATA/tenCate2002velocityRe1.5.csv", delimiter=",", names=["X", "Y"])


#ax[0, 0].plot(data_3_RondonDense['X'], data_3_RondonDense['Y'], linestyle='none', marker='o',color="k",fillstyle="none")



#version=input("choose version:")
version="../RestingSphereOverset/sphereAndBackground3D/"
sol_Mixt_RefSedim_A = version

time_t=0

datatime_Mixt = np.genfromtxt(sol_Mixt_RefSedim_A+'time',delimiter='\t', names=True)
dataVz_Mixt = np.genfromtxt(sol_Mixt_RefSedim_A+'vz',delimiter='\t', names=True)
dataPosz_Mixt = np.genfromtxt(sol_Mixt_RefSedim_A+'zcenter',delimiter='\t', names=True)


mixtUz_RefSedim_A=[]
mixtT_RefSedim_A=[]
mixtPosz_RefSedim_A=[]

count=0
for i in range(10,len(dataVz_Mixt)-10):
	if i%1==0:
		mixtUz_RefSedim_A.append(dataVz_Mixt[i][0])
		mixtT_RefSedim_A.append(datatime_Mixt[i][0]-time_t)
		mixtPosz_RefSedim_A.append(dataPosz_Mixt[i][0]/H)






#########################################################################


alhpaV=0.9
	
import matplotlib.pyplot as plt




################################################
limX=1


fig, axs = plt.subplots(2)
axs[0].set_title('Sphere position')
axs[0].plot(mixtT_RefSedim_A, mixtPosz_RefSedim_A,color='orange',alpha=alhpaV,label="SedFoam")
axs[0].plot(data_TenCate_displacement['X'], data_TenCate_displacement['Y'], linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[0].set_xlim([-0, limX])
axs[0].set_ylim([ -1.51,0])
axs[0].set_xlabel('t [s]')
axs[0].set_ylabel('z/D [-]')
axs[0].grid()

axs[1].set_title('Sphere velocity')
axs[1].plot(mixtT_RefSedim_A, mixtUz_RefSedim_A,color='orange',alpha=alhpaV,label="SedFoam")
axs[1].plot(data_TenCate_vel['X'], data_TenCate_vel['Y'], linestyle='none', marker='o',color="k",fillstyle="none",label="ten Cate (2002)")
axs[1].set_xlabel('t [s]')
axs[1].set_ylabel('Vertical velocity [m/s]')
axs[1].set_xlim([-0, limX])
axs[1].set_ylim([ -0.0401,0.0001])
axs[1].grid()
axs[1].legend(loc="upper right",fontsize="7")
fig.tight_layout()
fig.savefig('Figures/FallingSphere.png')
		
plt.show()



