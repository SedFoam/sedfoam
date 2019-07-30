#
# Import section
#
import subprocess
import sys
import numpy as np
import fluidfoam
from pylab import figure, subplot, axis, xlabel, ylabel, show, savefig, plot
from pylab import title, matplotlib
import matplotlib.gridspec as gridspec
import matplotlib as mpl
#
# Change fontsize
#
matplotlib.rcParams.update({'font.size': 20})
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10

#
# Change subplot sizes
#
gs = gridspec.GridSpec(2, 2)
gs.update(left=0.1, right=0.95, top=0.95,
          bottom=0.075, wspace=0.125, hspace=0.125)
#
# Figure size
#
figwidth = 12
figheight = 6
figheight2 = 12

#########################################
#
# Loading experimental results
#
#execfile('DATA/exp_lmsgc.py')
exec(open("DATA/exp_lmsgc.py").read())
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
    #
    # FIGURE 1: Interface positions Vs Time
    #
    figure(num=1, figsize=(figwidth, figheight),
           dpi=60, facecolor='w', edgecolor='w')

    plot(t_pvb + tadj, zint_pvb + 0.1, 'ob',
         t_pvb + tadj, zint2_pvb + 0.1, 'or')
    plot(time, zint, '-b', time, zint2, '-r')
    ylabel('y (m)')
    xlabel('t (s)')
    axis([0, tmax, zmin, zmax])

    show(block=False)
    savefig('Figures/res1_tuto1.png', facecolor='w',
            edgecolor='w', format='png')

#    toto=raw_input('hit a key to continue')
#
# Figure 2: Concentration profiles
#
figure(num=2, figsize=(figwidth, figheight2),
       dpi=60, facecolor='w', edgecolor='w')

for i in np.arange(4):
    if i == 0:
        ax = subplot(gs[0, 0])
    elif i == 1:
        ax = subplot(gs[0, 1])
    elif i == 2:
        ax = subplot(gs[1, 0])
    elif i == 3:
        ax = subplot(gs[1, 1])

    iexp = 7 * i + 1
    titi = np.where(time[:] >= t_pvb[0][iexp] + tadj)
    if np.size(titi) == 0:
        inum = Nt - 1
    else:
        inum = titi[0][0]
    print('texp= ' + str(t_pvb[0][iexp] + tadj) + 's - tnum='+ str(time[inum]) + 's')

    ax.plot(alphat[:, inum], Y[:], '-r',
            as_pvb[:, iexp], z_pvb[:, iexp] + 0.1, '--b')
    title('t=' + str(t_pvb[0][iexp] + tadj) + 's')
    axis([-0.05, 0.65, zmin, zmax])
    if (i == 0) or (i == 2):
        ylabel('y (m)')
    else:
        ax.set_yticklabels([''])

    if (i == 2) or (i == 3):
        xlabel(r'$\alpha$')
    else:
        ax.set_xticklabels([''])


savefig('Figures/res2_tuto1.png', facecolor='w', edgecolor='w', format='png')

show(block=False)
# Fix Python 2.x.
try: input = raw_input
except NameError: pass
tto = input("Hit a key to close the figure")
