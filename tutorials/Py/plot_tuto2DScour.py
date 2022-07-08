# =======================================================================
#                         GENERAL INFORMATION
# =======================================================================
""" Script plotting the evolution of the bed interface versus time and
    the 2D snapshots of the bed evolution at differents time.
"""
#
# ---------------- Module General Import and Declarations ---------------
#
import numpy as np
import fluidfoam
from pylab import matplotlib, xlabel, ylabel, axis, title, xticks, yticks, show
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#
# --------------------Figure definition----------------------
#
# Change fontsize
#
matplotlib.rcParams.update({'font.size': 12})
# Figure size
#
figwidth = 18
figheight = 9
#
gs1 = gridspec.GridSpec(1, 1)
gs1.update(left=0.08, right=0.9, top=0.95, bottom=0.1,
           wspace=0.175, hspace=0.2)
gs2 = gridspec.GridSpec(2, 2)
gs2.update(left=0.08, right=0.9, top=0.95, bottom=0.1,
           wspace=0.175, hspace=0.2)
#
# ---------------------Reading part---------------------
#


def readOpenFoam(sol, t0, Nt, Dt, Nx, Ny, Nz, N):
    tlist = t0 + np.arange(Nt) * Dt
    timeRange = [(repr(item).rstrip('0')).rstrip('.') for item in tlist]

    Xd, Yd, Zd = fluidfoam.readmesh(sol)  # ,shape=(Nx,Ny,Nz))
    toto = np.where(Xd >= 0)
    shape = (Nx, Ny, Nz)

    X = np.array((Nx, Ny, Nz))
    Y = np.array((Nx, Ny, Nz))
    # if the shape is prescribed, reshape the arrays
    if (max(shape) != 1):
        X = np.reshape(Xd[toto], shape, order="F")
        Y = np.reshape(Yd[toto], shape, order="F")
    Xp = X[:, 0, 0]
    Yp = Y[0, :, 0]
    time = np.zeros(Nt)
    alphap = np.zeros((Nx, Ny, Nz, Nt))
    ybed = np.zeros((Nx, Nt))  # bed interface

    k = -1
    for t in timeRange:
        print("Reading time: " + str(t) + " s")
        k = k + 1
        alphad = fluidfoam.readscalar(
            sol, t + '/', 'alpha_a')
        # if the shape is prescribed, reshape the arrays
        if (max(shape) != 1):
            alpha = np.reshape(alphad[toto], shape, order="F")
        alphap[:, :, :, k] = alpha[:, :, :]

        time[k] = float(t)
        for i in range(Nx):
            ybed[i, k] = Y[
                i, np.max((np.where(alphap[i, :, 0, k] > 0.57))), 0]

    return alphap, Xp, Yp, tlist, ybed

# Parameters
t0 = 0.
Nt = 13
Dt = 5


h0 = 0.15

#
# Call the reading function
#
print("-------------------------")
print("        case 2DScour     ")
print("-------------------------")
Nx = 1000
Ny = 264
Nz = 1
N = Nx * Ny
readu = 0
basepath = '../RAS/'
figpath = './Figures/'
figName = 'scour_evolTime_komega_MuI.png'
figName2 = 'scour_snapshots_komega_MuI.png'

case = '2DScour'

sol = basepath + case + '/'
alphap, Xp, Yp, tlist, ybed = readOpenFoam(sol, t0, Nt, Dt, Nx, Ny, Nz, N)

# -------------------------------------------------------------
#                          plot section
# --------------------------------------------------------------
#
# -----------fig1: evolution of the bed interface versus time--------------
#
fig = plt.figure(num=1, figsize=(18, 8), dpi=100, facecolor='w', edgecolor='w')
ax0 = fig.add_subplot(gs1[0, 0])
tplot = [0, 2, 4, 6]
nplot = -1
sym = ['--+', '-.', '--', '.']
for k in range(len(tplot)):
    nplot = nplot + 1
    l0 = ax0.plot(Xp, ybed[:, tplot[k]], 'k' + sym[nplot],
                  linewidth=2, label=r't=' + str(tplot[nplot] * Dt) + ' s')
ax0.legend(numpoints=1, loc='best', fontsize=20)
xlabel('x (m)')
ylabel('interface position (m)')
axis([0, 0.5, -0.02, 0.01])

plt.savefig(figpath+figName,format='png',dpi = 300)

#
# -----------fig2: 2D bed interface at several time--------------
#
fig = plt.figure(num=2, figsize=(16, 16),
                 dpi=100, facecolor='w', edgecolor='w')

nlevels = 100
levels = [0, 0.05, 0.1, 0.15, 0.2, 0.25,
          0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.625]
levels2 = [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

tfig = [0, 2]
tfig2 = [4, 6]
nplot = -1
for k in range(len(tfig)):
    nplot = nplot + 1
    alphat = np.transpose(alphap[:, :, 0, tfig[k]])
    alphat2 = np.transpose(alphap[:, :, 0, tfig2[k]])

    ax1 = fig.add_subplot(gs2[nplot, 0])
    title('t=' + str(tfig[nplot] * Dt) + ' s')
    CS_colors = plt.contourf(Xp, Yp, alphat, nlevels,  cmap=plt.cm.coolwarm)
    l0 = ax1.plot(Xp, ybed[:, tfig[k]], '--k', linewidth=1)
    axis([0, 0.5, -0.03, 0.03])
    xticks([], [])
    props = dict(boxstyle='round', facecolor='none', alpha=1.0)

    xticks([], [])
    if nplot == 1:
        xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5],
               ['0', '0.1', '0.2', '0.3', '0.4', '0.5'])
        xlabel('x (m)')

    yticks([-0.025, 0, 0.025], ['-0.025', '0', '0.025'])

    ax2 = fig.add_subplot(gs2[nplot, 1])
    title('t=' + str(tfig2[nplot] * Dt) + ' s')
    CS_colors = plt.contourf(Xp, Yp, alphat2, nlevels,  cmap=plt.cm.coolwarm)
    l0 = ax2.plot(Xp, ybed[:, tfig2[k]], '--k', linewidth=1)
    axis([0, 0.5, -0.03, 0.03])
    xticks([], [])
    yticks([], [])
    props = dict(boxstyle='round', facecolor='none', alpha=1.0)
    if nplot == 1:
        xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5],
               ['0', '0.1', '0.2', '0.3', '0.4', '0.5'])
        xlabel('x (m)')

    yticks([-0.025, 0, 0.025], ['-0.025', '0', '0.025'])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.925, 0.1, 0.02, 0.75])
cb = plt.colorbar(CS_colors, cax=cbar_ax)

labels = ['0', '0.2', '0.4', '0.6']
loc = [0, 0.2, 0.4, 0.6]
cb.set_ticks(loc)
cb.set_ticklabels(labels)

cb.set_label(r'$\alpha$')


plt.savefig(figpath+figName2,format='png',dpi = 300)

show()
