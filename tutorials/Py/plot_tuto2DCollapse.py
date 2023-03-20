import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import math
import os
from struct import unpack
import subprocess
from fluidfoam import readvector, readscalar, readmesh, readprobes
import csv
import sys
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rc('font', family='serif')

#########################################
# Extracting data from Rondon's experiments
#########################################

# Deposit morphology

Rondon3D = np.genfromtxt("DATA/Rondon2012/RondonDenseL6H4.2_T3s.csv", delimiter="\t", names=["X", "Y"])
Rondon6D = np.genfromtxt("DATA/Rondon2012/RondonDenseL6H4.2_T6s.csv", delimiter="\t", names=["X", "Y"])
Ron12D = np.genfromtxt("DATA/Rondon2012/RondonDenseL6H4.2_T12s.csv", delimiter="\t", names=["X", "Y"])

RondonAL = np.genfromtxt("DATA/Rondon2012/RondonLooseL6H4.8_T0.66s.csv", delimiter="\t", names=["X", "Y"])
RondonBL = np.genfromtxt("DATA/Rondon2012/RondonLooseL6H4.8_T1.32s.csv", delimiter="\t", names=["X", "Y"])
RondonCL = np.genfromtxt("DATA/Rondon2012/RondonLooseL6H4.8_T4s.csv", delimiter="\t", names=["X", "Y"])

# Evolution of  pore pressure with time

data1 = np.genfromtxt('DATA/Rondon2012/RondonDenseL6H4.2_FRONT', delimiter='\t', names=True)
time_exp_frontDense = data1['time']
front_L_expDense = data1['length']

data1 = np.genfromtxt('DATA/Rondon2012/RondonLooseL6H4.8_FRONT', delimiter='\t', names=True)
time_exp_frontLoose = data1['time']
front_L_expLoose = data1['length']

data1 = np.genfromtxt('DATA/Rondon2012/RondonDenseL6H4.2_Pressure', delimiter='\t', names=True)
time_exp_pressureDense = data1['time']
pressure_expDense = data1['pressure']

data1 = np.genfromtxt('DATA/Rondon2012/RondonLooseL6H4.8_Pressure', delimiter='\t', names=True)
time_exp_pressureLoose = data1['time']
pressure_expLoose = data1['pressure']

#########################################
# Data From SedFoam
#########################################

# Plots the contour of sediment concentration
levels = np.arange(-0.01, 0.45, 0.4)

# LOOSE CASE

TimeList = ["0.66", "1.32", "3.96"]
sol_loose = '../laminar/2DCollapse/loose/'

isExist = os.path.exists(sol_loose+"/postProcessing")
if not isExist:
    print("It seems you have not run the initially loose granular collapse :( ")
    sys.exit()

probeLoose, time_loose, PorePressure_loose = readprobes(sol_loose, time_name='mergeTime', name='p_rbgh')

x_l, y_l, z_l = readmesh(sol_loose, structured=True)
nx_l, ny_l, nz_l = x_l.shape

PorePressure_loose_2cm = []
PorePressure_loose_3cm = []
for i in range(len(time_loose)):
    PorePressure_loose_2cm.append(PorePressure_loose[i, 0][0])
    PorePressure_loose_3cm.append(PorePressure_loose[i, 1][0])

alpha_loose = []
for i in range(len(TimeList)):
    alpha_loose.append(readscalar(sol_loose, TimeList[i], 'alpha.a', structured=True))

# DENSE CASE

TimeList = ["3", "6", "12"]
sol_dense = '../laminar/2DCollapse/dense/'

isExist = os.path.exists(sol_dense+"/postProcessing")
if not isExist:
    print("It seems you have not run the initially dense granular collapse :( ")
    sys.exit()

probeDense, time_dense, PorePressure_dense = readprobes(sol_dense, time_name='mergeTime', name='p_rbgh')

x_d, y_d, z_d = readmesh(sol_dense, structured=True)
nx_d, ny_d, nz_d = x_d.shape

time_dense_v_k1_A_list = []
PorePressure_dense_2cm = []
PorePressure_dense_3cm = []
for i in range(len(time_dense)):
    PorePressure_dense_2cm.append(PorePressure_dense[i, 0][0])
    PorePressure_dense_3cm.append(PorePressure_dense[i, 1][0])

alpha_dense=[]
for i in range(len(TimeList)):
    alpha_dense.append(readscalar(sol_dense, TimeList[i], 'alpha.a', structured=True))

# Morphology plot

fig, ax = plt.subplots(3, 2, figsize=(12, 6))
plt.subplots_adjust(hspace=None)
ax[0, 0].text(0.075, 0.058, 'Dense', fontsize=14, rotation=0)
ax[0, 1].text(0.08, 0.058, 'Loose', fontsize=14, rotation=0)

ax[0, 0].contour(x_d[:, :, nz_d//2], y_d[:, :, nz_d//2], alpha_dense[0][:, :, nz_d//2], colors='b', levels=levels)
ax[1, 0].contour(x_d[:, :, nz_d//2], y_d[:, :, nz_d//2], alpha_dense[1][:, :, nz_d//2], colors='b', levels=levels)
ax[2, 0].contour(x_d[:, :, nz_d//2], y_d[:, :, nz_d//2], alpha_dense[2][:, :, nz_d//2], colors='b', levels=levels)

ax[0, 0].plot(Rondon3D['X'], Rondon3D['Y'], linestyle='none', marker='o', c='k', fillstyle="none")
ax[1, 0].plot(Rondon6D['X'], Rondon6D['Y'], linestyle='none', marker='o', c='k', fillstyle="none")
ax[2, 0].plot(Ron12D['X'], Ron12D['Y'], linestyle='none', marker='o', c='k', fillstyle="none", label="Rondon (2012)")

lim_x = 0.20
lim_y = 0.055

ax[0, 0].set_title('t = 3s', x=0.15, y=1, pad=-14, fontsize=12)
ax[0, 0].set_xlim([0, lim_x])
ax[0, 0].set_ylim([0, lim_y])
ax[0, 0].set_ylabel('y [m]')
ax[0, 0].label_outer()
ax[0, 0].set_aspect('equal')

ax[1, 0].set_title('t = 6s', x=0.15, y=1, pad=-14, fontsize=12)
ax[1, 0].set_xlim([0, lim_x])
ax[1, 0].set_ylim([0, lim_y])
ax[1, 0].set_ylabel('y [m]')
ax[1, 0].label_outer()
ax[1, 0].set_aspect('equal')

ax[2, 0].set_title('t = 12s', x=0.15, y=1, pad=-14, fontsize=12)
ax[2, 0].set_xlim([0, lim_x])
ax[2, 0].set_ylim([0, lim_y])
ax[2, 0].set_xlabel('x [m]')
ax[2, 0].label_outer()
ax[2, 0].set_ylabel('y [m]')
ax[2, 0].set_aspect('equal')

ax[0, 1].contour(x_l[:, :, nz_l//2], y_l[:, :, nz_l//2], alpha_loose[0][:, :, nz_l//2], colors='r', levels=levels)
ax[1, 1].contour(x_l[:, :, nz_l//2], y_l[:, :, nz_l//2], alpha_loose[1][:, :, nz_l//2], colors='r', levels=levels)
ax[2, 1].contour(x_l[:, :, nz_l//2], y_l[:, :, nz_l//2], alpha_loose[2][:, :, nz_l//2], colors='r', levels=levels)

ax[0, 1].plot(RondonAL['X'], RondonAL['Y'], linestyle='none', marker='o', c='k', fillstyle="none")
ax[1, 1].plot(RondonBL['X'], RondonBL['Y'], linestyle='none', marker='o', c='k', fillstyle="none")
ax[2, 1].plot(RondonCL['X'], RondonCL['Y'], linestyle='none', marker='o', c='k', fillstyle="none")

ax[0, 1].set_title('t = 0.66s', x=0.15, y=1, pad=-14, fontsize=12)
ax[0, 1].set_xlim([0, lim_x])
ax[0, 1].set_ylim([0, lim_y])
ax[0, 1].set_ylabel('y [m]')
ax[0, 1].label_outer()
ax[0, 1].set_aspect('equal')

ax[1, 1].set_title('t = 1.32s', x=0.15, y=1, pad=-14, fontsize=12)
ax[1, 1].set_xlim([0, lim_x])
ax[1, 1].set_ylim([0, lim_y])
ax[1, 1].set_ylabel('y [m]')
ax[1, 1].label_outer()
ax[1, 1].set_aspect('equal')

ax[2, 1].set_title('t = 4s', x=0.15, y=1, pad=-14, fontsize=12)
ax[2, 1].set_xlim([0, lim_x])
ax[2, 1].set_ylim([0, lim_y])
ax[2, 1].set_xlabel('x [m]')
ax[2, 1].set_ylabel('y [m]')
ax[2, 1].label_outer()
ax[2, 1].set_aspect('equal')

plt.subplots_adjust(hspace=0.01, wspace=0.13)

plt.savefig('Figures/tuto2DCollapse_Morphology.png')

# PorePressure plot

plt.figure(figsize=(8, 3.5))
ExpDense, = plt.plot(time_exp_pressureDense, pressure_expDense, marker='o', fillstyle='full', linewidth = 0, color='k')
ExpLoose, = plt.plot(time_exp_pressureLoose, pressure_expLoose, marker='o', fillstyle='none', linewidth = 0, color='k')

plt.plot(time_loose, PorePressure_loose[:, 0], linestyle='-', color='r')
plt.plot(time_loose, PorePressure_loose[:, 1], linestyle='--', color='r')
shade_loose = plt.fill_between(time_loose, PorePressure_loose_2cm, PorePressure_loose_3cm, facecolor='r', alpha=0.5)

plt.plot(time_dense, PorePressure_dense[:, 0], linestyle='-', color='b', label='Dense - 2cm')
plt.plot(time_dense, PorePressure_dense[:, 1], linestyle='--', color='b', label='Dense- 3cm')
shade_dense = plt.fill_between(time_dense, PorePressure_dense_2cm, PorePressure_dense_3cm, facecolor='b', alpha=0.5)

plt.legend([ExpLoose, ExpDense, shade_loose, shade_dense], ["Exp. - Loose", "Exp. - Dense", "Loose", "Dense"], loc=0)

plt.xlabel('Time [$s$]', fontsize=15)
plt.ylabel('Excess of pore pressure [$Pa$]', fontsize=15)
plt.xlim(0, 20.01)
plt.ylim(-160, 150)
plt.tight_layout()
plt.savefig('Figures/tuto2DCollapse_PorePressure.png')

plt.show()
