
# Import section
#
import sys
import os
import scipy.io
import numpy as np
import fluidfoam
from pylab import *
from scipy.io.matlab import mio
import matplotlib.gridspec as gridspec
#
def rms(x):
    return np.sqrt(x.dot(x)/x.size)
#
# Change fontsize
#
matplotlib.rcParams.update({'font.size': 15})
# Figure size
#
figwidth=18
figheight=9



#
#
# 
#---------------Loading OpenFoam results--------------------
#
basepath='../'


#
# Loading OpenFoam results
#
casedir='1DBoundaryLayer/'
tout='1000'

sol=basepath+casedir

x, z, y = fluidfoam.readmesh(sol)
k      = fluidfoam.readscalar(sol, tout,'k')
U      = fluidfoam.readvector(sol,tout,'Ub')
Tauf   = fluidfoam.readtensor(sol, tout, 'Taub')
u=U[0,:]
#########################################
#
# Physical parameters
#

rhof=1
nu=7.2727e-5

wallShear = np.max(Tauf[3,:])/rhof

H=np.max(z)
Umax=np.max(U)
Um=np.trapz(u,z)/H

print(' Reb=',Um*H/nu,' Um=',Um,' m/s' )
utau=np.sqrt(np.max(np.abs(wallShear)))
print(' Re*=',utau*H/nu,' u*=',utau,' m/s' )
# 
#----------Loading literature results-------------------
#
#
# Original k-w Wilcox model 
#
npzfiles= np.load('DATA/kw_wilcox.npz')

uw         = npzfiles['arr_0']
zuw        = npzfiles['arr_1']
kw         = npzfiles['arr_2']
zkw        = npzfiles['arr_3']
Ruww       = npzfiles['arr_4']
zRuww      = npzfiles['arr_5']
epsilonw   = npzfiles['arr_6']
zpepsilonw = npzfiles['arr_7']
#prodw      = npzfiles['arr_8']
#zpprodw    = npzfiles['arr_9']

#
# Low-Reynolds k-w Wilcox model 
#
npzfiles= np.load('DATA/kwLowRe_wilcox.npz')

uwlr         = npzfiles['arr_0']
zuwlr        = npzfiles['arr_1']
kwlr         = npzfiles['arr_2']
zkwlr        = npzfiles['arr_3']
Ruwwlr       = npzfiles['arr_4']
zRuwwlr      = npzfiles['arr_5']
epsilonwlr   = npzfiles['arr_6']
zpepsilonwlr = npzfiles['arr_7']
#prodwlr      = npzfiles['arr_8']
#zpprodwlr    = npzfiles['arr_9']

#
# Original k-w Wilcox model 
#
npzfiles= np.load('DATA/DNS_mansour.npz')

uDNS         = npzfiles['arr_0']
zuDNS        = npzfiles['arr_1']
kDNS         = npzfiles['arr_2']
zkDNS        = npzfiles['arr_3']
RuwDNS       = npzfiles['arr_4']
zRuwDNS      = npzfiles['arr_5']
epsilonDNS   = npzfiles['arr_6']
zpepsilonDNS = npzfiles['arr_7']
#prodDNS      = npzfiles['arr_8']
#zpprodDNS    = npzfiles['arr_9']

#
#---------------------Figures----------------------
#


#############fig3
gs3 = gridspec.GridSpec(1,3)
fig = plt.figure(num=3,figsize=(18,8),dpi=60, facecolor='w', edgecolor='w')
ax = fig.add_subplot(gs3[0,0])
pO   = ax.plot(u/Umax,z/H,'-or',label="OpenFOAM")
pw   = ax.plot(uw,zuw,'xb',label="Wilcox (1998)")
pwlr = ax.plot(uwlr,zuwlr,'+g',label="Wilcox Low-Re (2006)")
pDNS = ax.plot(uDNS,zuDNS,'^k',label="DNS (1999)")
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles,labels,numpoints=1,loc='best')
xlabel(r'$u/u_{max}$',fontsize=20)
ylabel(r'z / h',fontsize=20)
ax.set_yscale('log')
ax.axis([0, 1.05, 1e-3, 1.02])

ax = fig.add_subplot(gs3[0,1])
pO   = ax.plot(u/Umax,z/H,'-or',label=r"OpenFOAM")
pw   = ax.plot(uw,zuw,'xb',label=r"Wilcox (1998)")
pwlr = ax.plot(uwlr,zuwlr,'+g',label=r"Wilcox Low-Re (2006)")
pDNS = ax.plot(uDNS,zuDNS,'^k',label=r"DNS (1999)")
xlabel(r'$u/u_{max}$',fontsize=20)
#ylabel(r'z / h')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels,numpoints=1,loc='best',fontsize=17)
ax.axis([0, 1.05, 0, 1.02])

ax = fig.add_subplot(gs3[0,2])
pO   = ax.plot(k/utau**2,z/H,'-or',label="OpenFOAM")
pw   = ax.plot(kw,zkw,'xb',label="Wilcox (1998)")
pwlr = ax.plot(kwlr,zkwlr,'+g',label="Wilcox Low-Re (20??)")
pDNS = ax.plot(kDNS,zkDNS,'^k',label="DNS (19??)")
xlabel(r'$k/u_{\tau}^2$',fontsize=20)
#ylabel(r'z / h')
ax.axis([0, 5, 0, 1.02])
show()



#u_interp = np.interp(zuw, z[:]/H,U[0,:]/Umax);
#rms_u = rms(u_interp - uw)
#assert(rms_u<=0.03)

#k_interp = np.interp(zkw,z[:]/H,k[:]/utau**2);
#rms_k = rms(k_interp - kw)
#assert(rms_k<=0.35)