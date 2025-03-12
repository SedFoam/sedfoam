#!usr/bin/env python3
import numpy as np
from pylab import figure, subplot, grid, axis, xlabel, ylabel, show
from pylab import xticks, savefig
import matplotlib.gridspec as gridspec
from matplotlib import rc
from pylab import matplotlib

import fluidfoam

font = {'family': 'serif', 'weight': 'normal', 'size': 18}
rc("text", usetex=True)
rc('font', **font)
matplotlib.rcParams.update({"font.size": 20})
#
gs = gridspec.GridSpec(2, 2)
gs.update(left=0.065, right=0.975, top=0.95, bottom=0.075, wspace=0.2, hspace=0.25)

# markersize
ms = 3

###############################################################################
# parameters
###############################################################################
h = 1
ReTau_dns = 180
nu = 3.6e-4
ustar_dns = ReTau_dns*nu/h
zmax = h*ustar_dns/nu

Umax = 20
Rmax = 1.2

# Load DNS chan180 data
chan_data = np.loadtxt('./DATA/DNS_channel_Retau_180.d', comments='##')

# Load LES 128 data obtained from OpenFoam
LES_data = np.loadtxt('./DATA/LES_channel_Retau_180.d')

case = "3DChannel180"
basepath = "../LES/"

sol = basepath + case + "/"
#
# Reading SedFoam results
#
tread = "120"

dpdx = 0.00395967
ustar = np.sqrt(dpdx/h)
print('ustar_dns=', ustar_dns, ' ustar_les=', ustar)

print("########## Reading openfoam data file ##########")
# Read vertical coordinates
x, y, z = fluidfoam.readmesh(sol, structured=True, precision=12)
ny = len(y[0, :, 0])
uny = int(ny / 2)
yi = y[0, 0:uny, 0] * ustar/nu

# Read temporaly averaged variables
ubf_ta = fluidfoam.readvector(sol, tread, "UbMeanF", structured=True, precision=12)
ubprimf_ta = fluidfoam.readtensor(sol, tread, "UbPrime2MeanF", structured=True, precision=12)
TKEMeanProd_ta = fluidfoam.readscalar(sol, tread, "TKEMeanProd_b", structured=True, precision=12)
SGSDissMean_ta = fluidfoam.readscalar(sol, tread, "SGSDissMean_b", structured=True, precision=12)
viscDissMean_ta = fluidfoam.readscalar(sol, tread, "viscDissMean_b", structured=True, precision=12)
turbDiffusionMean_ta = fluidfoam.readscalar(sol, tread, "turbDiffusionMean_b", structured=True, precision=12)

# symm tensor :   XX , XY , XZ , YY , YZ , ZZ
#                  0    1    2    3    4    5
# Averaging
ubf_ta = ubf_ta[:, :, 0:uny, :]
ubprimf_ta = ubprimf_ta[:, :, 0:uny, :]
TKEMeanProd_ta = TKEMeanProd_ta[:, 0:uny, :]
SGSDissMean_ta = SGSDissMean_ta[:, 0:uny, :]
viscDissMean_ta = viscDissMean_ta[:, 0:uny, :]
turbDiffusionMean_ta = turbDiffusionMean_ta[:, 0:uny, :]

# Spatial averaging of temporaly averaged variables
ubf_a = np.mean(np.mean(ubf_ta, 3), 1)/ustar
ubprimf_a = np.mean(np.mean(ubprimf_ta, 3), 1)/ustar**2
TKE_a = 0.5*(ubprimf_a[0, :]+ubprimf_a[4, :]+ubprimf_a[8, :])
TKEMeanProd_a = np.mean(np.mean(TKEMeanProd_ta, 2), 0)*nu/ustar**4
SGSDissMean_a = np.mean(np.mean(SGSDissMean_ta, 2), 0)*nu/ustar**4
viscDissMean_a = np.mean(np.mean(viscDissMean_ta, 2), 0)*nu/ustar**4
turbDiffusionMean_a = np.mean(np.mean(turbDiffusionMean_ta, 2), 0)*nu/ustar**4

###############################################################################
# Plotting
###############################################################################
fig = figure(num=1, figsize=(12, 6), dpi=100, facecolor="w", edgecolor="w")

# Velocity profile
ax00 = subplot(gs[0, 0])
ax00.plot(chan_data[:, 4], chan_data[:, 0], 'ok', markersize=ms,
          label=r'$DNS$')
ax00.plot(LES_data[:, 4], LES_data[:, 0], '--b', lw=1,
          label=r'$LES (64^3)$ ref')
ax00.plot(ubf_a[0, :], yi[:], '-b', lw=2,
          label=r'$LES (32^3)$')

ax00.set_ylabel("$y^+$")
ax00.set_xlabel("$U/u_*$")
ax00.set_xlim(0, Umax)
ax00.set_ylim(1e-1, zmax)
ax00.set_yscale("log")
ax00.grid()
ax00.legend(loc='lower right', ncol=1, fontsize=14)

###############################################################################
ax01 = subplot(gs[0, 1])
ax01.plot(-chan_data[:, 5], chan_data[:, 0], 'ok', markersize=ms,
          label=r'$DNS$')
ax01.plot(-LES_data[:, 5], LES_data[:, 0], '--b', lw=1,
          label=r'$LES (128^3)$ ref')
ax01.plot(-ubprimf_a[1, :], yi[:], '-b', lw=2,
          label=r"$LES (32^3)$")
ax01.set_ylabel("$y^+$")
ax01.set_xlabel("$R_{xy}/u_*^2$")
ax01.set_xlim(0, 1)
ax01.set_ylim(0, zmax)
ax01.grid()
ax01.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax01.legend(fontsize=14)

###############################################################################
ax10 = subplot(gs[1, 0])
ax10.plot(chan_data[:, 6], chan_data[:, 0], 'ok', markersize=ms,
          label=r'$DNS$')
ax10.plot(LES_data[:, 6], LES_data[:, 0], '--b', lw=1,
          label=r'$LES (128^3)$ ref')
ax10.plot(TKE_a[:], yi[:], '-b', lw=2, label=r"$LES  (32^3)$")
ax10.set_ylabel("$y^+$")
ax10.set_xlabel("$TKE/u_*^2$")
ax10.set_xlim(0, 4.5)
ax10.set_ylim(0, zmax)
ax10.grid()
ax10.set_xticks([0, 1.5, 3, 4.5])
ax10.legend(fontsize=14)

###############################################################################
ax11 = subplot(gs[1, 1])
ax11.plot(chan_data[:, 1], chan_data[:, 0], 'ok', markersize=ms,
          label=r'$DNS[\epsilon]$')
ax11.plot(LES_data[:, 1], LES_data[:, 0], '--k', lw=1,
          label=r'$LES[\epsilon]$ ref')
ax11.plot(SGSDissMean_a + viscDissMean_a, yi[:], '-k', lw=2,
          label=r"$LES[\epsilon_{Res} + \epsilon_{Sgs}]  (32^3)$")
ax11.plot(viscDissMean_a, yi[:], '-m', lw=1,
          label=r"$LES[\epsilon_{Res}]  (32^3)$")

ax11.plot(chan_data[:, 3], chan_data[:, 0], 'or', markersize=ms,
          label=r'$DNS[TurbDiff]$')
ax11.plot(LES_data[:, 3], LES_data[:, 0], '--r', lw=1,
          label=r'$LES[TurbDiff]$ ref')
ax11.plot(turbDiffusionMean_a, yi[:], '-r', lw=2,
          label=r'$LES[TurbDiff] (32^3)$')

ax11.plot(chan_data[:, 2], chan_data[:, 0], 'ob', markersize=ms,
          label=r'$DNS[Prod]$')
ax11.plot(LES_data[:, 2], LES_data[:, 0], '--b',  lw=1,
          label=r'$LES[Prod]$ ref')
ax11.plot(TKEMeanProd_a, yi[:], '-b', lw=2,
          label=r'$LES[Prod]  (32^3)$')

ax11.set_ylabel("$y^+$")
ax11.set_ylim(0, zmax)
ax11.set_xlim(-0.25, 0.25)
ax11.grid()
ax11.legend(fontsize=8)
fig.text(0.70, 0.015, 'Loss', ha='center')
fig.text(0.85, 0.015, 'Gain', ha='center')

savefig('Fig/tuto3DChannel180_FavreAverage.png', format='png')

show(block=True)
