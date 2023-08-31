import numpy as np
import subprocess
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from fluidfoam import MeshVisu,readmesh,readscalar
import csv 
import os

degrees=180
L=0.15
H=0.12
nx=150
ny=120
dx=L/nx
dy=H/ny
dz=L*degrees*np.pi/180
VolAproxCell=dx*dy*dz*180
density=2500
time_down=0.15
alphaSusp=5e-3

############# SedFoam #################
sol = "../2DAxisymmetricMovingPlate/Upwards/"
try:
    proc = subprocess.Popen(
        ["foamListTimes", "-latestTime", "-case", sol], stdout=subprocess.PIPE
    )
except FileNotFoundError:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
tread = output.decode().rstrip().split("\n")[0]

myMesh = MeshVisu( path =  sol)
centr,vol=myMesh.getVolume()
Sus_material_laminar=[]
time=[]
steps=0.01
final=float(tread)
timeLine_laminar=np.arange(steps,final+steps,steps)
for i in range(len(timeLine_laminar)):
	tot_mat=0
	count=1
	timeStr=str(round(timeLine_laminar[i],2))
	if os.path.exists(sol+timeStr)==1:			
		alpha = readscalar(sol, timeStr, 'alpha.a')
		for j in range(len(alpha)):
			if (alpha[j]<alphaSusp):
				tot_mat+=alpha[j]*vol[j]*density*1000*degrees
		Sus_material_laminar.append(tot_mat)
		time.append(float(timeLine_laminar[i]/time_down))

plt.figure()
plt.plot(time,Sus_material_laminar,linestyle='--', linewidth=2, color='b',label='overSedDymFoam')
plt.legend(prop={'size':7},loc=1)
plt.xlabel('Time/$T_{lag}$ [$-$]', fontsize=15)
plt.ylabel('Bead weight [$g$]', fontsize=15)
plt.tight_layout()
plt.savefig('Figures/plotSuspensionVStime_motvingDiskUpwards.png', dpi=200)
plt.show()
