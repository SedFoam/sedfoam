
#!usr/bin/env python3
from netCDF4 import Dataset
import fluidfoam
import numpy as np
import os
import numpy as np
import fluidfoam
from fluidfoam import OpenFoamSimu
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import*
from diagnostic_function import*


A = 10
E = 3.5
F = 3

###########

mpl.style.use('classic')

plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.size'] =3
plt.rcParams['ytick.minor.width'] = 1

rc('axes',linewidth=1)

font = {'family' : 'serif', 'weight' : 'normal', 'size' : A}
plt.rc('font', **font)
plt.rc('text', usetex=True)

###########


plt.ion()

##########################################
# parameters
#
#############################
h = 1


case = 180#None
ReTau_dns = 180
nu = 3.6e-4

ustar_dns = ReTau_dns*nu/h

zmax = h*ustar_dns/nu


Umax = 25
Rmax = 1.2

################### DNS chan180 data

# load chan180 data

chan_data = np.loadtxt('DNS_channel_Retau_180.d', comments = '##')

###############   Load LES 64 data obtained from OpenFoam 

LES_data = np.loadtxt('LES_channel_Retau_180.d')


#################################
### Plotting
###########################
fig, ax = plt.subplots(2,2, figsize = (2*E,2*F))


ax[0][0].plot(chan_data[:, 4], chan_data[:, 0],                 'k', lw = 1, label = r'$Chan$', dashes = (4,2))

ax[0][0].plot(LES_data[:, 4], LES_data[:, 0],    'g', lw = 1, label = r'$LES$')

ax[0][0].set_ylabel("$y^+$")
ax[0][0].set_xlabel("$U/u_*$")
ax[0][0].set_xlim(0,Umax)
ax[0][0].set_ylim(1e-1,zmax)
ax[0][0].set_yscale("log")
ax[0][0].grid()
#ax[0][0].legend(bbox_to_anchor= (1.0, 1.0), ncol=1, fontsize = 0.7*A)
ax[0][0].legend(loc = 'lower right', ncol=1, fontsize = 0.7*A)
fig.tight_layout()

#######################################

ax[0][1].plot(-chan_data[:, 5], chan_data[:, 0],     'k', lw = 1, label = r'$Chan$', dashes = (4,2))

ax[0][1].plot(-LES_data[:, 5], LES_data[:, 0],  'g', lw = 1, label = r"$LES$")

ax[0][1].set_ylabel("$y^+$")
ax[0][1].set_xlabel("$R_{xy}/u_*^2$")
ax[0][1].set_xlim(0,0.8)
ax[0][1].set_ylim(0,zmax)
ax[0][1].grid()
ax[0][1].set_xticks([0, 0.2, 0.4, 0.8])
ax[0][1].legend(bbox_to_anchor= (1.0, 1.0), ncol=1, fontsize = 0.7*A)
fig.tight_layout()

#################################################
ax[1][0].plot(chan_data[:, 6], chan_data[:, 0],     'k', lw = 1, label = r'$Chan$', dashes = (4,2))

ax[1][0].plot(LES_data[:, 6], LES_data[:, 0], 'g', lw = 1, label = r"$LES$")	

ax[1][0].set_ylabel("$y^+$")
ax[1][0].set_xlabel("$TKE/u_*^2$")
ax[1][0].set_xlim(0, 4.5)
ax[1][0].set_ylim(0,zmax)
ax[1][0].grid()
ax[1][0].set_xticks([0, 1.5, 3, 4.5])
ax[1][0].legend(bbox_to_anchor= (1.0, 1.0), ncol=1, fontsize = 0.7*A)
fig.tight_layout()

#################

#############################

ax[1][1].plot(chan_data[:, 1], chan_data[:, 0],     'k', lw = 1, label = r'$Chan[\epsilon]$', dashes = (4,2))
ax[1][1].plot(LES_data[:, 1], LES_data[:, 0], 'g', lw = 1, label = r"$sum[Res + Sgs]$")
ax[1][1].plot(chan_data[:, 3] , chan_data[:, 0],                 'r', markersize = 3, lw = 1, label = r'$Chan[TurbDiff]$', dashes = (4,2))
ax[1][1].plot(LES_data[:, 3],  LES_data[:, 0],    '*r', markersize = 3,  lw = 1, label = r'$LES[TurbDiff]$')
ax[1][1].plot(chan_data[:, 2] , chan_data[:, 0],                 'b', markersize = 4, lw = 1, label = r'$Chan[Prod]$', dashes = (4,2))
ax[1][1].plot(LES_data[:, 2], LES_data[:, 0],    '+b', markersize = 4,  lw = 1, label = r'$LES[Prod]$' )

ax[1][1].set_ylabel("$y^+$")
#ax[1][1].set_xticks([-0.18, -0.12, -0.06,0.0])
ax[1][1].set_ylim(0,zmax)
#ax.set_ylabel("$P/\epsilon$")
ax[1][1].set_xlim(-0.25, 0.25)
#ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
ax[1][1].grid()
ax[1][1].legend(bbox_to_anchor= (1.0, 1.0), ncol=1, fontsize = 0.55*A)
fig.text(0.70, 0.015, 'Loss', ha='center')
fig.text(0.85, 0.015, 'Gain', ha='center')

fig.tight_layout()


plt.savefig('3DChannel180_FavreAverage.pdf')
plt.savefig('3DChannel180_FavreAverage.png')

plt.show()
