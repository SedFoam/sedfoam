from netCDF4 import Dataset
import subprocess
import fluidfoam
import numpy as np
import sys


case = '3DChannel560'
basepath = '../'
sol = basepath + case + '/'

#
# Reading SedFoam results
#
tread = '8.5'

print('########## Writing averaged data file ##########')
# Read vertical coordinates
x, y, z = fluidfoam.readmesh(sol, True, precision=12)
ny = len(y[0, :, 0])
uny = int(ny/2)
yi = y[0, 0:uny, 0]

# Read temporaly averaged variables
alpha_ta = fluidfoam.readscalar(sol, tread, 'alpha_aMean', True, precision=12)
ubf_ta = fluidfoam.readvector(sol, tread, 'UbMeanF', True, precision=12)
uaf_ta = fluidfoam.readvector(sol, tread, 'UaMeanF', True, precision=12)
ubprimf_ta = fluidfoam.readtensor(sol, tread, 'UbPrime2MeanF', True,
                                  precision=12)
uaprimf_ta = fluidfoam.readtensor(sol, tread, 'UaPrime2MeanF', True,
                                  precision=12)

# Usable data
alpha_ta = alpha_ta[:, 0:uny, :]
ubf_ta = ubf_ta[:, :, 0:uny, :]
uaf_ta = uaf_ta[:, :, 0:uny, :]
ubprimf_ta = ubprimf_ta[:, :, 0:uny, :]
uaprimf_ta = uaprimf_ta[:, :, 0:uny, :]

# Spatial averaging of temporaly averaged variables
alpha_a = np.mean(np.mean(alpha_ta, 2), 0)
ubf_a = np.mean(np.mean(ubf_ta, 3), 1)
uaf_a = np.mean(np.mean(uaf_ta, 3), 1)
ubprimf_a = np.mean(np.mean(ubprimf_ta, 3), 1)
uaprimf_a = np.mean(np.mean(uaprimf_ta, 3), 1)

# NetCDF file creation
rootgrp = Dataset(sol+'/postProcessing/'+case+'_averaged.nc', 'w')

# Dimensions creation
rootgrp.createDimension('vect', 3)
rootgrp.createDimension('tens', 9)
rootgrp.createDimension('coord', uny)

# Variables creation
y_file = rootgrp.createVariable('y', np.float64, 'coord')
alpha_a_file = rootgrp.createVariable('alpha_a', np.float64, 'coord')
ubprimf_a_file = rootgrp.createVariable('ubprimf_a', np.float64, ('tens',
                                                                  'coord'))
uaprimf_a_file = rootgrp.createVariable('uaprimf_a', np.float64, ('tens',
                                                                  'coord'))
ubf_a_file = rootgrp.createVariable('ubf_a', np.float64, ('vect', 'coord'))
uaf_a_file = rootgrp.createVariable('uaf_a', np.float64, ('vect', 'coord'))

# Writing variables
y_file[:] = yi
alpha_a_file[:] = alpha_a
ubprimf_a_file[:, :] = ubprimf_a
uaprimf_a_file[:, :] = uaprimf_a
ubf_a_file[:, :] = ubf_a
uaf_a_file[:, :] = uaf_a

# File closing
rootgrp.close()