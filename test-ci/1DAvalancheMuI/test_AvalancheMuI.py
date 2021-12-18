import analyticBagnold
import subprocess
import os
import numpy as np
import fluidfoam
import sys

def rms(x):
    return np.sqrt(x.dot(x)/x.size)

#
# compute the analytical solution in dimensionless form
#
# dimensional parameters
d = 0.02
rho_f = 1e-3
rho_p = 1
drho = rho_p - rho_f
g = 1

# dimensionless parameters
mus = 0.38
mu2 = 0.64
I0 = 0.3
Bphi = 0.31
phi0 = 0.59999999

beta = 26*np.pi/180.


#########################################
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
    print("Did you loaded OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
tread = output.decode().rstrip().split('\n')[0]

#########################################
# Reading SedFoam results
#########################################

X,Y,Z = fluidfoam.readmesh(sol)
alpha = fluidfoam.readscalar(sol, tread, 'alpha.a')
Ua = fluidfoam.readvector(sol, tread, 'U.a')
Ub = fluidfoam.readvector(sol, tread, 'U.b')
pff = fluidfoam.readscalar(sol, tread, 'pff')
pa = fluidfoam.readscalar(sol, tread, 'pa')
muI = fluidfoam.readscalar(sol, tread, 'muI')
nuEffa = fluidfoam.readscalar(sol, tread, 'nuEffa')
nuEffb = fluidfoam.readscalar(sol, tread, 'nuEffb')
nuFra = fluidfoam.readscalar(sol, tread, 'nuFra')
Tauf = fluidfoam.readtensor(sol, tread, 'Taub')
Taus = fluidfoam.readtensor(sol, tread, 'Taua')
try:
    gradUa = fluidfoam.readtensor(sol, tread, 'grad(U.a)')
except:
    print("grad(Ua) was not found -> Introduce - postProcess -func 'grad(U.a)' - in the command line")
    os.system("postProcess -func \'grad(Ua)\' -time "+tread)
    gradUa = fluidfoam.readtensor(sol, tread, 'grad(U.a)')


Ny = np.size(Y)
U = np.zeros(Ny)
U = alpha[:] * Ua[0, :] + (1 - alpha[:]) * Ub[0, :]

taua = np.zeros(Ny)
taub = np.zeros(Ny)
taua = Taus[3, :]
taub =  Tauf[3, :]
if (np.size(pff))==1:
    pff=pff*np.ones(Ny)

print("max(Ua)=" + str(np.amax(Ua)) + " m/s")

dudy = gradUa[3,:]
nu = np.zeros(Ny)
for i in range(Ny - 1):
    nu[i] = taua[0]/(dudy[i]*rho_p)

granularBed = np.where(pa+pff<(rho_p-rho_f)*g*d/10)
nx = np.size(granularBed)
H=Y[np.min(granularBed)]/d

print("H/d=",H)
zmax = np.max(Y)/d

[xex,I,alphaex,uex,pex,muIex,tauex,duexdz,nuex] = \
             analyticBagnold.analyticBagnold(nx,H,g,d,rho_p,rho_f,phi0,I0,Bphi,mus,mu2,beta)
print("max(uex)=" + str(np.max(uex)) + " m/s" + "max(phi)=",np.max(alphaex),"max(nuex)=",np.max(nuex))


# =============================================================================

alpha_interp = np.interp(xex, Y[:],alpha);
rms_alpha = rms(alpha_interp - alphaex)
assert(rms_alpha<=0.05)

u_interp = np.interp(xex, Y[:],Ua[0,:]);
rms_u = rms(u_interp - uex)
assert(rms_u<=0.09)

muI_interp = np.interp(xex,Y[:], muI);
rms_muI = rms(muI_interp - muIex)
assert(rms_muI<=0.02)



