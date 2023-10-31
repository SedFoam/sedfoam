from fluidfoam.readpostpro import readforce
import numpy as np
import fluidfoam
import os
import sys

rhoSolid=1141
rhof=1000
R=0.0075
V=4*np.pi*R**3/3
g=-9.81
weigth=V*g*rhoSolid
weigthBouyancy=V*g*rhof
valuesChose=100
print ("Volume: ",V )
print ("bouyancy contr: ",weigthBouyancy )
print ("weight contr: ",weigth )
print ("submerged weight: ",(weigth -weigthBouyancy))
print ("reduced mass: ",(weigth -weigthBouyancy)/9.81)

############## loading overSedDymFoam results ###########
sol="../RestingSphereMorphing/"
if os.path.exists(sol+'time')==False:
        print ("To generate the postprocessing files you need to execute makeFiles in the case first")   
        sys.exit()
datatime_Mixt_A = np.genfromtxt(sol+'time.txt',delimiter='\t', names=True)
dataVz_Mixt_A = np.genfromtxt(sol+'vz.txt',delimiter='\t', names=True)
dataPosz_Mixt_A = np.genfromtxt(sol+'zcenter.txt',delimiter='\t', names=True)
data_ForceFluidz_Mixt_A = np.genfromtxt(sol+'pressureFluid_z.txt',delimiter='\t', names=True)
data_ForceViscz_Mixt_A = np.genfromtxt(sol+'pressureVisc_z.txt',delimiter='\t', names=True)
data_ForceParticlez_Mixt_A = np.genfromtxt(sol+'pressureParticle_z.txt',delimiter='\t', names=True)

mixtUz_RefSedim_A=[]
mixtT_RefSedim_A=[]
mixtPosz_RefSedim_A=[]
mixtFfluid_z_RefSedim_A=[]
mixtFVisc_z_RefSedim_A=[]
mixtFPart_z_RefSedim_A=[]
mixtTOTAL_z_RefSedim_A=[]

count=0
for i in range(10,len(dataVz_Mixt_A)-10):
	if i%valuesChose==0:
		mixtUz_RefSedim_A.append(dataVz_Mixt_A[i][0])
		mixtT_RefSedim_A.append(datatime_Mixt_A[i][0])
		mixtPosz_RefSedim_A.append( dataPosz_Mixt_A[i][0]/(2*R))
		mixtFfluid_z_RefSedim_A.append(data_ForceFluidz_Mixt_A[i*2+1][0])
		mixtFVisc_z_RefSedim_A.append(data_ForceViscz_Mixt_A[i*2+1][0])
		mixtFPart_z_RefSedim_A.append(data_ForceParticlez_Mixt_A[i*2+1][0])
		mixtTOTAL_z_RefSedim_A.append(data_ForceFluidz_Mixt_A[i*2+1][0]+data_ForceParticlez_Mixt_A[i*2+1][0]+data_ForceViscz_Mixt_A[i*2+1][0]+weigth-weigthBouyancy)

#########################################################################

alhpaV=0.9
	
from pylab import matplotlib, plt, show

plt.rcParams.update({'font.size': 16})
plt.rc('font', family='serif')
plt.rc('text',usetex=True)
font = {'family':'serif','size':16, 'serif': ['computer modern roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':16})
matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']
matplotlib.rcParams["legend.framealpha"] = None

AWTime=[0,5]
AW=[weigth-weigthBouyancy,weigth-weigthBouyancy]

plt.rcParams["figure.figsize"] = (6,6)

fig, axs = plt.subplots(3, gridspec_kw=dict(height_ratios=[0.5, 0.5,1]))
#axs[0].set_title('Sphere position')
axs[0].plot(mixtT_RefSedim_A, mixtPosz_RefSedim_A,color='royalblue', linestyle='-')
axs[0].set_xlim([-0, mixtT_RefSedim_A[-1]])
axs[0].set_xticklabels([])
axs[0].set_ylabel('z/D [-]',fontsize="16")
axs[0].grid()

#axs[1].set_title('Sphere velocity')
axs[1].plot(mixtT_RefSedim_A, mixtUz_RefSedim_A,color='royalblue', linestyle='-')
axs[1].set_xticklabels([])
axs[1].set_ylabel('$v_{sphere}$ [m/s]',fontsize="16")
axs[1].set_xlim([-0, mixtT_RefSedim_A[-1]])
axs[1].grid()

#axs[2].set_title('Total force')
axs[2].plot(mixtT_RefSedim_A, mixtFPart_z_RefSedim_A,color='sienna',linestyle="--",alpha=alhpaV,label="Particle pressure contribution")
axs[2].plot(mixtT_RefSedim_A, mixtFfluid_z_RefSedim_A,color='g',linestyle="-.",alpha=alhpaV,label="Fluid pressure contribution")
axs[2].plot(mixtT_RefSedim_A, mixtFVisc_z_RefSedim_A,color='b',linestyle="--",alpha=alhpaV,label="Viscous contribution")
axs[2].plot(AWTime, AW,color='gray',alpha=alhpaV,linestyle="--",label=" Apparent weight")
axs[2].plot(mixtT_RefSedim_A, mixtTOTAL_z_RefSedim_A,color='k',alpha=alhpaV,label="Total force")
axs[2].set_xlim([-0, mixtT_RefSedim_A[-1]])
axs[2].set_xlabel('Time [s]')
axs[2].set_ylabel('Force [N]')
axs[2].grid()
axs[2].legend(loc="upper left",fontsize="12" ,framealpha=0.5)

fig.tight_layout()
fig.savefig('Figures/RestingSphereDynamicMeshComparison.png', dpi=300)
	
plt.show()


