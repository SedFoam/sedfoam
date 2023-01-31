from netCDF4 import Dataset
import fluidfoam
import numpy as np
import os


case = "3DOscillSheetFlow"
basepath = "../LES/"
sol = basepath + case + "/"

#
# Reading SedFoam results
#

print("########## Writing averaged data file ##########")
x, y, z = fluidfoam.readmesh(sol, True, precision=12)
yi = y[0, :, 0]

dir_list = os.listdir(sol)
time_list = []

for directory in dir_list:
    try:
        float(directory)
        time_list.append(directory)
    except:
        pass
time_list.sort(key=float)
time_list = np.array(time_list)

write_alpha = True
write_ub = True
write_ua = True

# Wave period
period = 5

# Number of divisions per period
n_div = 20

#
# Read data over four periods and average
#

# alpha_a
if write_alpha:
    rootgrp = Dataset(sol + "/postProcessing/" + case + "_alpha.nc", "w")
    rootgrp.createDimension("n", len(yi))
    rootgrp.createDimension("p", n_div)
    pos_file = rootgrp.createVariable("pos", np.float64, "n")
    phase_file = rootgrp.createVariable("phase", np.float64, "p")
    alpha_file = rootgrp.createVariable("alpha", np.float64, ("n", "p"))
    pos_file[:] = yi
    phase_file[:] = np.linspace(0, 5, n_div)
    for i, time in enumerate(time_list[1:n_div]):
        alpha1 = fluidfoam.readscalar(sol, time, "alpha_a", True, precision=12)
        alpha2 = fluidfoam.readscalar(
            sol, time_list[i + n_div + 1], "alpha_a", True, precision=12
        )
        alpha3 = fluidfoam.readscalar(
            sol, time_list[i + 2 * n_div + 1], "alpha_a", True, precision=12
        )
        alpha4 = fluidfoam.readscalar(
            sol, time_list[i + 3 * n_div + 1], "alpha_a", True, precision=12
        )
        alpha_a = (
            np.mean(np.mean(alpha1, 2), 0)
            + np.mean(np.mean(alpha2, 2), 0)
            + np.mean(np.mean(alpha3, 2), 0)
            + np.mean(np.mean(alpha4, 2), 0)
        ) / 4
        alpha_file[:, i] = alpha_a
    rootgrp.close()

# Ub / UbPrim
if write_ub:
    rootgrp_u = Dataset(sol + "/postProcessing/" + case + "_Ub.nc", "w")
    rootgrp_uu = Dataset(sol + "/postProcessing/" + case + "_UbPrim.nc", "w")
    rootgrp_u.createDimension("n", len(yi))
    rootgrp_u.createDimension("p", n_div)
    rootgrp_uu.createDimension("n", len(yi))
    rootgrp_uu.createDimension("p", n_div)
    rootgrp_u.createDimension("vect", 3)
    rootgrp_uu.createDimension("prim", 4)
    pos_file_u = rootgrp_u.createVariable("pos", np.float64, "n")
    pos_file_uu = rootgrp_uu.createVariable("pos", np.float64, "n")
    phase_file_u = rootgrp_u.createVariable("phase", np.float64, "p")
    phase_file_uu = rootgrp_uu.createVariable("phase", np.float64, "p")
    ub_file = rootgrp_u.createVariable("Ub", np.float64, ("vect", "n", "p"))
    ubprim_file = rootgrp_uu.createVariable("UbPrim", np.float64, ("prim", "n", "p"))
    pos_file_u[:] = yi
    phase_file_u[:] = np.linspace(0, 5, n_div)
    pos_file_uu[:] = yi
    phase_file_uu[:] = np.linspace(0, 5, n_div)
    for i, time in enumerate(time_list[1:n_div]):
        ub1 = fluidfoam.readvector(sol, time, "Ub", True, precision=12)
        ub2 = fluidfoam.readvector(
            sol, time_list[i + n_div + 1], "Ub", True, precision=12
        )
        ub3 = fluidfoam.readvector(
            sol, time_list[i + 2 * n_div + 1], "Ub", True, precision=12
        )
        ub4 = fluidfoam.readvector(
            sol, time_list[i + 3 * n_div + 1], "Ub", True, precision=12
        )
        ub_a = (
            np.mean(np.mean(ub1, 3), 1)
            + np.mean(np.mean(ub2, 3), 1)
            + np.mean(np.mean(ub3, 3), 1)
            + np.mean(np.mean(ub4, 3), 1)
        ) / 4
        ub_file[:, :, i] = ub_a
        ub_abox = np.zeros(ub1.shape)
        for j, xi in enumerate(x[:, 0, 0]):
            for k, zi in enumerate(z[0, 0, :]):
                ub_abox[:, j, :, k] = ub_a
        ubprim1 = ub1 - ub_abox
        ubprim2 = ub2 - ub_abox
        ubprim3 = ub3 - ub_abox
        ubprim4 = ub4 - ub_abox
        urms_a = np.sqrt(
            np.mean(
                np.mean(
                    (
                        ubprim1[0, :, :, :] * ubprim1[0, :, :, :]
                        + ubprim2[0, :, :, :] * ubprim2[0, :, :, :]
                        + ubprim3[0, :, :, :] * ubprim3[0, :, :, :]
                        + ubprim4[0, :, :, :] * ubprim4[0, :, :, :]
                    )
                    / 4,
                    2,
                ),
                0,
            )
        )
        vrms_a = np.sqrt(
            np.mean(
                np.mean(
                    (
                        ubprim1[1, :, :, :] * ubprim1[1, :, :, :]
                        + ubprim2[1, :, :, :] * ubprim2[1, :, :, :]
                        + ubprim3[1, :, :, :] * ubprim3[1, :, :, :]
                        + ubprim4[1, :, :, :] * ubprim4[1, :, :, :]
                    )
                    / 4,
                    2,
                ),
                0,
            )
        )
        wrms_a = np.sqrt(
            np.mean(
                np.mean(
                    (
                        ubprim1[2, :, :, :] * ubprim1[2, :, :, :]
                        + ubprim2[2, :, :, :] * ubprim2[2, :, :, :]
                        + ubprim3[2, :, :, :] * ubprim3[2, :, :, :]
                        + ubprim4[2, :, :, :] * ubprim4[2, :, :, :]
                    )
                    / 4,
                    2,
                ),
                0,
            )
        )
        uv_a = np.mean(
            np.mean(
                (
                    ubprim1[0, :, :, :] * ubprim1[1, :, :, :]
                    + ubprim2[0, :, :, :] * ubprim2[1, :, :, :]
                    + ubprim3[0, :, :, :] * ubprim3[1, :, :, :]
                    + ubprim4[0, :, :, :] * ubprim4[1, :, :, :]
                )
                / 4,
                2,
            ),
            0,
        )
        ubprim_file[0, :, i] = urms_a
        ubprim_file[1, :, i] = vrms_a
        ubprim_file[2, :, i] = wrms_a
        ubprim_file[3, :, i] = uv_a
    rootgrp_u.close()
    rootgrp_uu.close()

# Ua / UaPrim
if write_ua:
    rootgrp_u = Dataset(sol + "/postProcessing/" + case + "_Ua.nc", "w")
    rootgrp_uu = Dataset(sol + "/postProcessing/" + case + "_UaPrim.nc", "w")
    rootgrp_u.createDimension("n", len(yi))
    rootgrp_u.createDimension("p", n_div)
    rootgrp_uu.createDimension("n", len(yi))
    rootgrp_uu.createDimension("p", n_div)
    rootgrp_u.createDimension("vect", 3)
    rootgrp_uu.createDimension("prim", 4)
    pos_file_u = rootgrp_u.createVariable("pos", np.float64, "n")
    pos_file_uu = rootgrp_uu.createVariable("pos", np.float64, "n")
    phase_file_u = rootgrp_u.createVariable("phase", np.float64, "p")
    phase_file_uu = rootgrp_uu.createVariable("phase", np.float64, "p")
    ua_file = rootgrp_u.createVariable("Ua", np.float64, ("vect", "n", "p"))
    uaprim_file = rootgrp_uu.createVariable("UaPrim", np.float64, ("prim", "n", "p"))
    pos_file_u[:] = yi
    phase_file_u[:] = np.linspace(0, 5, n_div)
    pos_file_uu[:] = yi
    phase_file_uu[:] = np.linspace(0, 5, n_div)
    for i, time in enumerate(time_list[1:n_div]):
        ua1 = fluidfoam.readvector(sol, time, "Ua", True, precision=12)
        ua2 = fluidfoam.readvector(
            sol, time_list[i + n_div + 1], "Ua", True, precision=12
        )
        ua3 = fluidfoam.readvector(
            sol, time_list[i + 2 * n_div + 1], "Ua", True, precision=12
        )
        ua4 = fluidfoam.readvector(
            sol, time_list[i + 3 * n_div + 1], "Ua", True, precision=12
        )
        ua_a = (
            np.mean(np.mean(ua1, 3), 1)
            + np.mean(np.mean(ua2, 3), 1)
            + np.mean(np.mean(ua3, 3), 1)
            + np.mean(np.mean(ua4, 3), 1)
        ) / 4
        ua_file[:, :, i] = ua_a
        ua_abox = np.zeros(ua1.shape)
        for j, xi in enumerate(x[:, 0, 0]):
            for k, zi in enumerate(z[0, 0, :]):
                ua_abox[:, j, :, k] = ua_a
        uaprim1 = ua1 - ua_abox
        uaprim2 = ua2 - ua_abox
        uaprim3 = ua3 - ua_abox
        uaprim4 = ua4 - ua_abox
        urms_a = np.sqrt(
            np.mean(
                np.mean(
                    (
                        uaprim1[0, :, :, :] * uaprim1[0, :, :, :]
                        + uaprim2[0, :, :, :] * uaprim2[0, :, :, :]
                        + uaprim3[0, :, :, :] * uaprim3[0, :, :, :]
                        + uaprim4[0, :, :, :] * uaprim4[0, :, :, :]
                    )
                    / 4,
                    2,
                ),
                0,
            )
        )
        vrms_a = np.sqrt(
            np.mean(
                np.mean(
                    (
                        uaprim1[1, :, :, :] * uaprim1[1, :, :, :]
                        + uaprim2[1, :, :, :] * uaprim2[1, :, :, :]
                        + uaprim3[1, :, :, :] * uaprim3[1, :, :, :]
                        + uaprim4[1, :, :, :] * uaprim4[1, :, :, :]
                    )
                    / 4,
                    2,
                ),
                0,
            )
        )
        wrms_a = np.sqrt(
            np.mean(
                np.mean(
                    (
                        uaprim1[2, :, :, :] * uaprim1[2, :, :, :]
                        + uaprim2[2, :, :, :] * uaprim2[2, :, :, :]
                        + uaprim3[2, :, :, :] * uaprim3[2, :, :, :]
                        + uaprim4[2, :, :, :] * uaprim4[2, :, :, :]
                    )
                    / 4,
                    2,
                ),
                0,
            )
        )
        uv_a = np.mean(
            np.mean(
                (
                    uaprim1[0, :, :, :] * uaprim1[1, :, :, :]
                    + uaprim2[0, :, :, :] * uaprim2[1, :, :, :]
                    + uaprim3[0, :, :, :] * uaprim3[1, :, :, :]
                    + uaprim4[0, :, :, :] * uaprim4[1, :, :, :]
                )
                / 4,
                2,
            ),
            0,
        )
        uaprim_file[0, :, i] = urms_a
        uaprim_file[1, :, i] = vrms_a
        uaprim_file[2, :, i] = wrms_a
        uaprim_file[3, :, i] = uv_a
    rootgrp_u.close()
    rootgrp_uu.close()
