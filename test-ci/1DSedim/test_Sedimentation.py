#
# Import section
#
import subprocess
import sys
import numpy as np
import fluidfoam

#########################################
#
# Loading experimental results
#
#execfile('DATA/exp_lmsgc.py')
exec(open("../../tutorials/Py/DATA/exp_lmsgc.py").read())
#########################################
# Loading OpenFoam results
#
#
#
#
case = '1DSedim'
basepath = '../'
sol = basepath + case + '/'

Nx = 1
Ny = 120
Nz = 1

eps_file = sol + case + '.eps'

#
# Reading SedFoam results
#
try:
    proc = subprocess.Popen(
        ['foamListTimes', '-case', sol], stdout=subprocess.PIPE)
except:
    print("foamListTimes : command not found")
    print("Do you have load OpenFoam environement?")
    sys.exit(0)
output = proc.stdout.read()
#tread = list(output.split('\n'))
tread = output.decode().rstrip().split('\n')


del tread[-1]
Nt = len(tread)
time = np.zeros(Nt)
X, Y, Z = fluidfoam.readmesh(sol)
alphat = np.zeros((Ny, Nt))

k = -1
for t in tread:
    print("Reading time: %s s" % t)
    k = k + 1

    alphat[:, k] = fluidfoam.readscalar(sol, t + '/', 'alpha_a')
    time[k] = float(t)


#
# parameter
#
zmin = 0.
zmax = np.max(Y)

tmax = 1800.
tadj = 172.

fontsize = 18.
#
# calcul zint et zint2
#
if Nt > 1:

    asint2 = 0.55
    asint = 0.25

    zint = np.zeros(Nt)
    zint2 = np.zeros(Nt)

    for i in np.arange(Nt):
        # zint
        toto = np.where(alphat[:, i] < asint)
        if np.size(toto) == 0:
            zint[i] = Y[Ny - 1]
        else:
            zint[i] = Y[toto[0][0]]
    # zint2
        toto2 = np.where(alphat[:, i] <= asint2)
        if np.size(toto2) == 0:
            zint2[i] = Y[0]
        else:
            zint2[i] = Y[toto2[0][0]]
 
    zint_pvb_interp = np.interp(t_pvb+tadj, time, zint);
    zint2_pvb_interp = np.interp(t_pvb+tadj, time, zint2);
    assert(np.allclose(zint_pvb +0.1, zint_pvb_interp, atol=1e-2))
    assert(np.allclose(zint2_pvb +0.1, zint2_pvb_interp, atol=1e-2))

