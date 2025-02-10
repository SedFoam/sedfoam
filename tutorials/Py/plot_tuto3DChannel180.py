
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
ReTau_les64 = 180
dpdx_les64  = 0.00391054 #0.00402415

ustar_dns = ReTau_dns*nu/h

zmax = h*ustar_dns/nu

Umax = 25
Rmax = 1.2

ustarN_les64  = (dpdx_les64*h)**(0.5)

TKEnorm_les64  = nu/ustarN_les64**4

################### DNS chan180 data

# load chan590 data

yUdns,yplusUdns,Udns,dum,dum,dum,dum = np.loadtxt('/.fsnet/data/legi/calcul9/home/sharma6ma/useful/project/17LES_SHEET/manohar/Sheet_Flow/Channel_Flow/data/DNS/chan_data/chan%d/profiles/chan%d.means'%(case, case),unpack=True)

yRdns,yplusRdns,Ruudns,Rvvdns,Rwwdns,Ruvdns,dum,dum = np.loadtxt('/.fsnet/data/legi/calcul9/home/sharma6ma/useful/project/17LES_SHEET/manohar/Sheet_Flow/Channel_Flow/data/DNS/chan_data/chan%d/profiles/chan%d.reystress'%(case, case),unpack=True)

ydns,yplusdns,dissipdns,producdns,pstraindns,pdiffdns,tdiffdns,vdiffdns,bal = np.loadtxt('/.fsnet/data/legi/calcul9/home/sharma6ma/useful/project/17LES_SHEET/manohar/Sheet_Flow/Channel_Flow/data/DNS/chan_data/chan%d/balances/chan%d.kbal'%(case, case),unpack=True)

################
###############   Load LES 64 data obtained from OpenFoam
##########

path = '/.fsnet/data/legi/calcul9/home/sharma6ma/useful/project/17LES_SHEET/manohar/Sheet_Flow/Channel_Flow/Numerical_data/LES/grid_64_sedfoam_validation/TKE_one_file/start_0'
time_frame = '2261.49711998406337'

prec = 9

# Read vertical coordinates

x, y, z = fluidfoam.readmesh(path, time_name='constant', structured=True, precision=prec)
ny = len(y[0, :, 0])
uny = int(ny)
yi = y[0, 0:uny, 0]

ub       = fluidfoam.readvector(path, time_frame, 'UbMeanF_b', True, precision = prec)
UbPrime2MeanF = fluidfoam.readtensor(path, time_frame, 'UbPrime2MeanF_b', True, precision = prec)
viscDissMean         = fluidfoam.readscalar(path, time_frame, 'viscDissMeanJ_b', True, precision=prec)
SGSDissMean         = fluidfoam.readscalar(path, time_frame, 'SGSDissMeanJ_b', True, precision=prec)
TKEMeanProd         = fluidfoam.readscalar(path, time_frame, 'TKEMeanProd_b', True, precision=prec)
turbDiffusionMean   = fluidfoam.readscalar(path, time_frame, 'turbDiffusionMean_b', True, precision=prec)

yb = yi

Ubx = ub[0]
Ubx = np.mean(np.mean(Ubx,2),0)

UbPrime2 = np.mean(np.mean(UbPrime2MeanF,3),1)

Rxy_les64 = UbPrime2[1,:]

viscDissMean = np.mean(np.mean(viscDissMean,2),0)
SGSDissMean  = np.mean(np.mean(SGSDissMean,2),0)
TKEProdMean = np.mean(np.mean(TKEMeanProd,2),0)
turbDiffusionMean = np.mean(np.mean(turbDiffusionMean,2),0)

###############
## compute Yplus
###############

yplus  = yb*ustarN_les64/nu

###############
## TKE
#############

tke_les64   = 0.5*(UbPrime2[0,:] + UbPrime2[4,:] + UbPrime2[8,:])

tkedns = 0.5*(Ruudns + Rvvdns + Rwwdns)

#################################
### Plotting
###########################
fig, ax = plt.subplots(2,2, figsize = (2*E,2*F))

ax[0][0].plot(Udns, yplusUdns,                 'k', lw = 1, label = r'$Chan$', dashes = (4,2))

ax[0][0].plot(Ubx/ustarN_les64, yplus,    'g', lw = 1, label = r'$LES$')

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

ax[0][1].plot(-Ruvdns, yplusUdns,     'k', lw = 1, label = r'$Chan$', dashes = (4,2))

ax[0][1].plot(-Rxy_les64/ustarN_les64**2, yplus,  'g', lw = 1, label = r"$LES$")

ax[0][1].set_ylabel("$y^+$")
ax[0][1].set_xlabel("$R_{xy}/u_*^2$")
ax[0][1].set_xlim(0,0.8)
ax[0][1].set_ylim(0,zmax)
ax[0][1].grid()
ax[0][1].set_xticks([0, 0.2, 0.4, 0.8])
ax[0][1].legend(bbox_to_anchor= (1.0, 1.0), ncol=1, fontsize = 0.7*A)
fig.tight_layout()

#################################################
ax[1][0].plot(tkedns, yplusUdns,     'k', lw = 1, label = r'$Chan$', dashes = (4,2))

ax[1][0].plot(tke_les64/ustarN_les64**2, yplus, 'g', lw = 1, label = r"$LES$")

ax[1][0].set_ylabel("$y^+$")
ax[1][0].set_xlabel("$TKE/u_*^2$")
ax[1][0].set_xlim(0, 4.5)
ax[1][0].set_ylim(0,zmax)
ax[1][0].grid()
ax[1][0].set_xticks([0, 1.5, 3, 4.5])
ax[1][0].legend(bbox_to_anchor= (1.0, 1.0), ncol=1, fontsize = 0.7*A)
fig.tight_layout()

#############################

ax[1][1].plot(dissipdns, yplusdns,     'k', lw = 1, label = r'$Chan[\epsilon]$', dashes = (4,2))
ax[1][1].plot((viscDissMean + SGSDissMean)*TKEnorm_les64, yplus, 'g', lw = 1, label = r"$sum[Res + Sgs]$")
ax[1][1].plot(tdiffdns , yplusdns,                 'r', markersize = 3, lw = 1, label = r'$Chan[TurbDiff]$', dashes = (4,2))
ax[1][1].plot(turbDiffusionMean*TKEnorm_les64, yplus,    '*r', markersize = 3,  lw = 1, label = r'$LES[TurbDiff]$')
ax[1][1].plot(producdns , yplusdns,                 'b', markersize = 4, lw = 1, label = r'$Chan[Prod]$', dashes = (4,2))
ax[1][1].plot(TKEProdMean*TKEnorm_les64, yplus,    '+b', markersize = 4,  lw = 1, label = r'$LES[Prod]$' )

ax[1][1].set_ylabel("$y^+$")
ax[1][1].set_ylim(0,zmax)
ax[1][1].set_xlim(-0.25, 0.25)
ax[1][1].grid()
ax[1][1].legend(bbox_to_anchor= (1.0, 1.0), ncol=1, fontsize = 0.55*A)
fig.text(0.70, 0.015, 'Loss', ha='center')
fig.text(0.85, 0.015, 'Gain', ha='center')

fig.tight_layout()

plt.savefig('3DChannel180_FavreAverage.pdf')
plt.savefig('3DChannel180_FavreAverage.png')

plt.show()

