#
# Import section
#
import numpy as np
import fluidfoam


def rms(x):
    return np.sqrt(x.dot(x) / x.size)


#
#
# ---------------Loading OpenFoam results--------------------
#
basepath = "./"
#
# Loading OpenFoam results
#
casedir = "./"
tout = "2500"

sol = basepath + casedir

x, z, y = fluidfoam.readmesh(sol)
k = fluidfoam.readscalar(sol, tout, "k.b")
U = fluidfoam.readvector(sol, tout, "U.b")
Tauf = fluidfoam.readtensor(sol, tout, "Taub")
u = U[0, :]
#########################################
#
# Physical parameters
#
rhof = 1
nu = 7.2727e-5

wallShear = np.max(Tauf[3, :]) / rhof

H = np.max(z)
Umax = np.max(U)
Um = np.trapz(u, z) / H

print(" Reb=", Um * H / nu, " Um=", Um, " m/s")
utau = np.sqrt(np.max(np.abs(wallShear)))
print(" Re*=", utau * H / nu, " u*=", utau, " m/s")
#
# ----------Loading literature results-------------------
#
#
# Original k-w Wilcox model
#
npzfiles = np.load("DATA/kw_wilcox.npz")

uw = npzfiles["arr_0"]
zuw = npzfiles["arr_1"]
kw = npzfiles["arr_2"]
zkw = npzfiles["arr_3"]
Ruww = npzfiles["arr_4"]
zRuww = npzfiles["arr_5"]
epsilonw = npzfiles["arr_6"]
zpepsilonw = npzfiles["arr_7"]
# prodw      = npzfiles['arr_8']
# zpprodw    = npzfiles['arr_9']

#
# Low-Reynolds k-w Wilcox model
#
npzfiles = np.load("DATA/kwLowRe_wilcox.npz")

uwlr = npzfiles["arr_0"]
zuwlr = npzfiles["arr_1"]
kwlr = npzfiles["arr_2"]
zkwlr = npzfiles["arr_3"]
Ruwwlr = npzfiles["arr_4"]
zRuwwlr = npzfiles["arr_5"]
epsilonwlr = npzfiles["arr_6"]
zpepsilonwlr = npzfiles["arr_7"]
# prodwlr      = npzfiles['arr_8']
# zpprodwlr    = npzfiles['arr_9']

#
# Original k-w Wilcox model
#
npzfiles = np.load("DATA/DNS_mansour.npz")

uDNS = npzfiles["arr_0"]
zuDNS = npzfiles["arr_1"]
kDNS = npzfiles["arr_2"]
zkDNS = npzfiles["arr_3"]
RuwDNS = npzfiles["arr_4"]
zRuwDNS = npzfiles["arr_5"]
epsilonDNS = npzfiles["arr_6"]
zpepsilonDNS = npzfiles["arr_7"]
# prodDNS      = npzfiles['arr_8']
# zpprodDNS    = npzfiles['arr_9']

u_interp = np.interp(zuw, z[:] / H, U[0, :] / Umax)
rms_u = rms(u_interp - uw)
print("rms_u=", rms_u)
assert rms_u <= 0.05

k_interp = np.interp(zkw, z[:] / H, k[:] / utau**2)
rms_k = rms(k_interp - kw)
print("rms_k=", rms_k)
assert rms_k <= 0.16
