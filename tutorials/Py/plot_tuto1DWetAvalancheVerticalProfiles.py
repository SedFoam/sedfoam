import numpy as np
import fluidfoam
from pylab import plt, matplotlib

plt.rcParams.update({'font.size': 16})
plt.rc('font', family='serif')

plt.rc('text',usetex=True)
font = {'family':'serif','size':16, 'serif': ['computer modern roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':16})
matplotlib.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']
matplotlib.rcParams["legend.framealpha"] = None

##################################################################
#Parameters to retrieve dimensionless variables
##################################################################

d= 160e-6 # particle diameter in m
gravity = 9.81 # gravity in m/s2
rhoFluid=1041  # fluid density in kg/m3
rhoSolid=2500  # solid density in kg/m3
h =  0.0049 # initial granular height in m
theta=25 # plane slope
val_p=(rhoSolid-rhoFluid)*gravity*h # pressure at the bottom
timeAdim=(d/gravity)**0.5
velAdim=1000.*(gravity*d)**0.5
pressureAdim=rhoFluid*h*gravity

#########################################
# Loading SedFoam results
#########################################

sol =  '../laminar/1DWetAvalanche'
X,Y,Z = fluidfoam.readmesh(sol)
tolAlpha=0.55

# this part of the script analyzes the vertical profiles at time=0 (before tilting the plane)
alpha_0 = fluidfoam.readscalar(sol, '200', 'alpha.a')
pff_0 = fluidfoam.readscalar(sol, '200', 'pff')
pa_0 = fluidfoam.readscalar(sol, '200', 'pa')
p_rbgh_0 = fluidfoam.readscalar(sol, '200', 'p_rbgh')
Ua_0= fluidfoam.readvector(sol,  '200', 'U.a')

vel_0=[]
phi_0=[]
y_0=[]
p_excess_0=[]
p_c_0=[]


for k in range(len(alpha_0)):
		if (alpha_0[k]>tolAlpha):
			vel_0.append(Ua_0[0, k]*1000/velAdim)
			phi_0.append(alpha_0[k])
			y_0.append(Y[k]/h)
			p_c_0.append(0)
			p_excess_0.append((p_rbgh_0[k]/val_p))
		else:
			break

# this part of the script analyzes the vertical profiles at times>0 (after tilting the plane)
times=[10,20,60,200,400]# please select specific times where data will be reconstructed

velocityProfiles=[]
phiProfiles=[]
yProfiles=[]
particlePressureProfiles=[]
excessPressureProfiles=[]

newHeight=h
for i in range(len(times)):
	time_v=times[i]
	tread=str(times[i]+200)+'/'
	alpha_A = fluidfoam.readscalar(sol, tread, 'alpha.a')
	Ua_A = fluidfoam.readvector(sol, tread, 'U.a')
	pff_A = fluidfoam.readscalar(sol, tread, 'pff')
	pa_A = fluidfoam.readscalar(sol, tread, 'pa')
	p_rbgh_A = fluidfoam.readscalar(sol, tread, 'p_rbgh')
	vel_values=[]
	phi_values=[]
	y_values=[]
	p_particle_values=[]
	vel_values=[]
	p_excess_values=[]
	
	for k in range(len(alpha_A)):
		if (alpha_A[k]<tolAlpha):
			newHeight=Y[k]
			break

	for k in range(len(alpha_A)):
		if (alpha_A[k]>tolAlpha and Y[k]<h):
			vel_values.append(Ua_A[0, k]*1000/velAdim)
			phi_values.append(alpha_A[k])
			y_values.append(Y[k]/h)
			p_particle_values.append((pff_A[k]+pa_A[k]-alpha_A[k]*(rhoSolid-rhoFluid)*gravity*np.cos(theta*np.pi/180)*(newHeight-Y[k]))/val_p)
			p_excess_values.append((p_rbgh_A[k]/val_p))
		else:
			break	

	velocityProfiles.append(vel_values)
	phiProfiles.append(phi_values)
	yProfiles.append(y_values)
	particlePressureProfiles.append(p_particle_values)
	excessPressureProfiles.append(p_excess_values)
	
	
	
			
#########################################
# 				Plots
#########################################
	
from matplotlib.ticker import StrMethodFormatter

# velocity profile
plt.figure()
plt.plot(vel_0,y_0, marker='o', markersize=0,linestyle='--', linewidth=1.5, color='k',label='$t=0s$')
plt.plot(velocityProfiles[0],yProfiles[0], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='powderblue')
plt.plot(velocityProfiles[1],yProfiles[1], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='deepskyblue')
plt.plot(velocityProfiles[2],yProfiles[2], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='dodgerblue')
plt.plot(velocityProfiles[3],yProfiles[3], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='royalblue')
plt.plot(velocityProfiles[4],yProfiles[4], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='navy')
plt.ylabel('$y/h_o$ [$-$]', fontsize=18)
plt.xlabel('$v^s/\\sqrt{gd}$ [$-$]', fontsize=18)
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}')) #
plt.grid()
plt.tight_layout()
#plt.legend(prop={'size':10.0},loc=0)
plt.savefig('Figures/VelocityProfile1D_phi0592.png', dpi=200)
	
	
# volume fraction profile
plt.figure()
plt.plot(phi_0,y_0, marker='o', markersize=0,linestyle='--', linewidth=1.5, color='k',label='$t=0s$')
plt.plot(phiProfiles[0],yProfiles[0], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='powderblue')
plt.plot(phiProfiles[1],yProfiles[1], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='deepskyblue')
plt.plot(phiProfiles[2],yProfiles[2], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='dodgerblue')
plt.plot(phiProfiles[3],yProfiles[3], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='royalblue')
plt.plot(phiProfiles[4],yProfiles[4], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='navy')
plt.ylabel('$y/h_o$ [$-$]', fontsize=18)
plt.xlabel('$\\phi$ [$-$]', fontsize=18)
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.3f}'))
plt.grid()
plt.tight_layout()
#plt.legend(loc=0)
plt.savefig('Figures/PhiProfile1D_phi0592.png', dpi=200)
	


# excess of fluid pore pressure profile
plt.figure()
plt.plot(p_excess_0,y_0, marker='o', markersize=0,linestyle='--', linewidth=1.5, color='k',label='$t=0s$')
plt.plot(excessPressureProfiles[0],yProfiles[0], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='powderblue',label='$t=%s s$'%times[0])
plt.plot(excessPressureProfiles[1],yProfiles[1], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='deepskyblue',label='$t=%s s$'%times[1])
plt.plot(excessPressureProfiles[2],yProfiles[2], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='dodgerblue',label='$t=%s s$'%times[2])
plt.plot(excessPressureProfiles[3],yProfiles[3], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='royalblue',label='$t=%s s$'%times[3])
plt.plot(excessPressureProfiles[4],yProfiles[4], marker='o', markersize=0,linestyle='-', linewidth=1.5, color='navy',label='$t=%s s$'%times[4])
plt.ylabel('$y/h_o$ [$-$]', fontsize=18)
plt.xlabel('$\\frac{p^f}{(\\rho^s - \\rho^f)g h_o}$ [$-$]', fontsize=21)
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
plt.grid()
plt.legend(prop={'size':15.0},loc=0)
plt.tight_layout()
plt.savefig('Figures/FluidPressureProfile1D_phi0592.png', dpi=200)

plt.show()
