import subprocess
import sys
import numpy as np
import fluidfoam

def rms(x):
    return np.sqrt(x.dot(x)/x.size)

#
zmin = 0
zmax = 1
#
# compute the analytical solution in dimensionless form
#
# dimensional parameters
d = 3.175
rho_f = 1e-3
rho_p = 1.
drho = rho_p - rho_f
g = 0.0

[phiexp,yexpphi] =np.genfromtxt('DATA/simpleShearVescovi2014_phiY.dat',delimiter=';',comments='#', unpack=True)
[uexp,yexp] =np.genfromtxt('DATA/simpleShearVescovi2014_uY.dat',delimiter=';',comments='#', unpack=True)
[thetaexp,yexptheta] =np.genfromtxt('DATA/simpleShearVescovi2014_TY.dat',delimiter=';',comments='#', unpack=True)

######################################
# Loading OpenFoam results
#########################################
case = './'
basepath = './'
sol = basepath + case + '/'

#
# Reading SedFoam results
#
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
X,Y,Z = fluidfoam.readmesh(sol,True,precision=prec)
nx,ny,nz = X.shape

alpha = fluidfoam.readscalar(sol, tread, 'alpha_a',True,precision=prec)
Ua    = fluidfoam.readvector(sol, tread, 'Ua',True,precision=prec)
pa    = fluidfoam.readscalar(sol, tread, 'pa',True,precision=prec)
Theta = fluidfoam.readscalar(sol, tread, 'Theta',True,precision=prec)

Ny = np.size(Y)
H = np.max(np.max(Y))
U = 1

iprof=0

phi_interp = np.interp(Y[0,:,iprof]/H,yexpphi, phiexp);
rms_phi = rms(phi_interp - alpha[0,:,iprof])
assert(rms_phi<=0.025)

u_interp = np.interp(Y[0,:,iprof]/H,yexp, uexp);
rms_u = rms(u_interp - Ua[0,0,:,iprof])
assert(rms_u<=0.026)

theta_interp = np.interp(Y[0,:,iprof]/H,yexptheta, thetaexp);
rms_theta = rms(theta_interp - Theta[0,:,iprof])/np.mean(thetaexp)
assert(rms_theta<=0.36)
