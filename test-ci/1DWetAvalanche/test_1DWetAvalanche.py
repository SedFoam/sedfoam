import subprocess
import numpy as np
import fluidfoam
import sys

def rms(x):
    return np.sqrt(x.dot(x)/x.size)


######################################
# Loading reference results
#########################################

data1 = np.genfromtxt('DATA/DATAphi592Num.txt',delimiter=' ', names=True)
pressure_592_num = data1['pressure']
velocity_592_num = data1['velocity']
alpha_592_num = data1['alpha']
dila_592_num = data1['dilatancy_angle']
y_592_num = data1['Y_data']

######################################
# Loading OpenFoam results
#########################################
case = './'
basepath = './'
sol = basepath + case + '/'


try:
    proc = subprocess.Popen(
        ['foamListTimes', '-latestTime', '-case', sol], stdout=subprocess.PIPE)
except:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
tread = output.decode().rstrip().split('\n')[0]

#########################################
# Reading SedFoam results
#########################################
prec=9
X,Y,Z = fluidfoam.readmesh(sol)

alpha_0 = fluidfoam.readscalar(sol, tread, 'alpha_a')
Ua_0   = fluidfoam.readvector(sol, tread, 'Ua')
p_rbgh_0    = fluidfoam.readscalar(sol, tread, 'p_rbgh')
delta_0 = fluidfoam.readscalar(sol, tread, 'delta')

Y_list=[]
alpha_list=[]
vel_list=[]
p_rbgh_list=[]
delta_list=[]

tolAlpha=0.54
for k, alpha0k in enumerate(alpha_0):
    if (alpha0k<tolAlpha):
        break
    Y_list.append(Y[k])
    alpha_list.append(alpha0k)
    vel_list.append(Ua_0[0,k])
    p_rbgh_list.append(p_rbgh_0[k])
    delta_list.append(delta_0[k])

#file0 = open("DATAphi592.txt","w")
#file0.writelines(["Y_data ","alpha ","velocity ","pressure ","dilatancy_angle "+'\n'])
#tolAlpha=0.54
#for k in range(len(alpha_0)):
		#if (alpha_0[k]<tolAlpha):
				#break
		#file0.writelines([str(Y[k])+" ",str(alpha_0[k])+" ",str(Ua_0[0,k])+" ",str(p_rbgh_0[k])+" ",str(delta_0[k])+'\n'])

iprof=0

# =============================================================================
phi_interp = np.interp(Y_list,y_592_num, alpha_592_num);
rms_phi = rms(phi_interp - alpha_list)
assert(rms_phi<=0.025)

# =============================================================================
pressure_interp = np.interp(Y_list,y_592_num, pressure_592_num);
rms_pressure = rms(pressure_interp - p_rbgh_list)
assert(rms_pressure<=3.1)

# =============================================================================
delta_interp = np.interp(Y_list,y_592_num, dila_592_num);
rms_delta = rms(delta_interp - delta_list)
assert(rms_delta<=0.025)

# =============================================================================
vel_interp = np.interp(Y_list,y_592_num, velocity_592_num);
rms_vel = rms(vel_interp - vel_list)
assert(rms_vel<=0.025)

