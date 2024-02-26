import numpy as np

rhoSolid=1141
rhof=1000
R=0.0075
V=4*np.pi*R**3/3
g=-9.81
weigth=V*g*rhoSolid
weigthBouyancy=V*g*rhof
valuesChose=10
print ("Volume: ",V )
print ("bouyancy contr: ",weigthBouyancy )
print ("weight contr: ",weigth )
print ("submerged weight: ",(weigth -weigthBouyancy))
print ("reduced mass: ",(weigth -weigthBouyancy)/9.81)


def moving_average(data, window_size):
    """
    Smooths the curve represented by a list of values using a moving average.

    Args:
    - data: List of values representing the curve.
    - window_size: Size of the window for the moving average.

    Returns:
    - smoothed_data: List of smoothed values.
    """
    smoothed_data = []
    for i in range(len(data)):
        start = max(0, i - window_size // 2)
        end = min(len(data), i + window_size // 2 + 1)
        avg = sum(data[start:end]) / (end - start)
        smoothed_data.append(avg)
    return smoothed_data




############## loading overSedDymFoam results ###########
sol="../RestingSphereMorphing/"
datatime_Mixt_A = np.genfromtxt(sol+'time.txt',delimiter='\t', names=True)
dataVz_Mixt_A = np.genfromtxt(sol+'vy.txt',delimiter='\t', names=True)
dataPosz_Mixt_A = np.genfromtxt(sol+'ycenter.txt',delimiter='\t', names=True)
data_ForceFluidz_Mixt_A = np.genfromtxt(sol+'pressureFluid_y.txt',delimiter='\t', names=True)
data_ForceViscz_Mixt_A = np.genfromtxt(sol+'pressureVisc_y.txt',delimiter='\t', names=True)
data_ForceParticlez_Mixt_A = np.genfromtxt(sol+'pressureParticle_y.txt',delimiter='\t', names=True)

mixtUz_RefSedim_A=[]
mixtT_RefSedim_A=[]
mixtPosz_RefSedim_A=[]
mixtFfluid_z_RefSedim_A=[]
mixtFVisc_z_RefSedim_A=[]
mixtFPart_z_RefSedim_A=[]
mixtTOTAL_z_RefSedim_A=[]

count=0
for i in range(100,len(dataVz_Mixt_A)-100):
	if i%valuesChose==0:
		mixtUz_RefSedim_A.append(dataVz_Mixt_A[i][0])
		mixtT_RefSedim_A.append(datatime_Mixt_A[i][0])
		mixtPosz_RefSedim_A.append( dataPosz_Mixt_A[i][0]/(2*R))
		mixtFfluid_z_RefSedim_A.append(data_ForceFluidz_Mixt_A[i*2+1][0])
		mixtFVisc_z_RefSedim_A.append(data_ForceViscz_Mixt_A[i*2+1][0])
		mixtFPart_z_RefSedim_A.append(data_ForceParticlez_Mixt_A[i*2+1][0])
		mixtTOTAL_z_RefSedim_A.append(data_ForceFluidz_Mixt_A[i*2+1][0]+data_ForceParticlez_Mixt_A[i*2+1][0]+data_ForceViscz_Mixt_A[i*2+1][0]+weigth-weigthBouyancy)

# Window size for the moving average
window = 10

# Smooth the curve
mixtUz_RefSedim_Ave = moving_average(mixtUz_RefSedim_A, window)
mixtPosz_RefSedim_Ave = moving_average(mixtPosz_RefSedim_A, window)
mixtFfluid_z_RefSedim_Ave = moving_average(mixtFfluid_z_RefSedim_A, window)
mixtFVisc_z_RefSedim_Ave = moving_average(mixtFVisc_z_RefSedim_A, window)
mixtFPart_z_RefSedim_Ave = moving_average(mixtFPart_z_RefSedim_A, window)
mixtTOTAL_z_RefSedim_Ave = moving_average(mixtTOTAL_z_RefSedim_A, window)
#########################################################################

alhpaV=1

from pylab import matplotlib, plt, show

plt.rcParams.update({'font.size': 16})
plt.rc('font', family='serif')
plt.rc('text',usetex=True)
font = {'family':'serif','size':16, 'serif': ['computer modern roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':16})
matplotlib.rcParams["legend.framealpha"] = None

AWTime=[0,0.31]
AW=[weigth-weigthBouyancy,weigth-weigthBouyancy]

plt.rcParams["figure.figsize"] = (9,6)

fig, axs = plt.subplots(3, gridspec_kw=dict(height_ratios=[0.5, 0.5,1]))
#axs[0].set_title('Sphere position')
axs[0].plot(mixtT_RefSedim_A, mixtPosz_RefSedim_Ave,color='royalblue', linestyle='-',label="Morphing")
axs[0].set_xlim([-0, AWTime[-1]])
axs[0].set_xticklabels([])
axs[0].set_ylabel('y/D [-]',fontsize="16")
axs[0].grid()

#axs[1].set_title('Sphere velocity')
axs[1].plot(mixtT_RefSedim_A, mixtUz_RefSedim_Ave,color='royalblue', linestyle='-',label="SedFoam - $\\phi=0$ - Relax=0.4")
axs[1].set_xticklabels([])
axs[1].set_ylabel('$v_{sphere}$ [m/s]',fontsize="16")
axs[1].set_xlim([-0, AWTime[-1]])
axs[1].grid()

#axs[2].set_title('Total force')
axs[2].plot(mixtT_RefSedim_A, mixtFPart_z_RefSedim_Ave,color='sienna',linestyle="--",alpha=alhpaV,label="Particle pressure contribution")
axs[2].plot(mixtT_RefSedim_A, mixtFfluid_z_RefSedim_Ave,color='g',linestyle="-.",alpha=alhpaV,label="Excess of fluid pressure contribution")
axs[2].plot(mixtT_RefSedim_A, mixtFVisc_z_RefSedim_Ave,color='b',linestyle="--",alpha=alhpaV,label="Viscous and frictional contribution")
axs[2].plot(AWTime, AW,color='gray',alpha=alhpaV,linestyle="--",label=" Apparent weight")
axs[2].plot(mixtT_RefSedim_A, mixtTOTAL_z_RefSedim_Ave,color='k',alpha=alhpaV,label="Total force")
axs[2].set_xlim([-0, AWTime[-1]])
axs[2].set_xlabel('Time [s]')
axs[2].set_ylabel('Force [N]')
axs[2].grid()
axs[2].legend(loc="upper left",fontsize="12" ,framealpha=0.5)

fig.tight_layout()
fig.savefig('Figures/RestingSphere.png', dpi=500)

plt.show()
