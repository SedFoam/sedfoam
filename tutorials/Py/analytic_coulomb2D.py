# --------------------------------------------------------------
#
#     AUTHOR : Julien CHAUCHAT                 DATE : 09/01/2010
#
# --------------------------------------------------------------
#
#     Function that solves the analytical solution proposed by
#       Ouriemi et al. (2009) for the case of the bed-load
#       transport in laminar shearing flows based on the
#       two-phase equations in the mixed-fluid form.
#
# --------------------------------------------------------------


def analytic_coulomb2D(nx, x, dpdx, hp, mus, phi0, eta_e):
    import numpy as np
# --------------------------------------------------------------
#
#  Calculation of the depth of the flowing layer
#
# --------------------------------------------------------------
    ratio = 1e0 / eta_e
    dum = - (1e0 - hp) / ratio * \
        (1e0 - (1e0 - ratio * dpdx / (mus * phi0 + dpdx)) ** 0.5)
    hc = hp - dum
# --------------------------------------------------------------
#
#  Calculation of the velocity profile
#
# --------------------------------------------------------------
    uex = np.zeros(nx)

    for i in range(nx):
        if x[i] < hc:
            #  ZONE III : Fixed granular layer
            uex[i] = 0e0
        elif x[i] <= hp:
            #  ZONE II  : Granular layer in motion (constant porosity)
            uex[i] = ratio * (dpdx + mus * phi0) * 0.5e0 * (x[i] - hc) ** 2
        elif x[i] > hp:
            #  ZONE I   : Pure fluid layer
            uex[i] = -dpdx * 0.5e0 * (1e0 - x[i]) * (x[i] - hp) \
                + 0.5e0 * ratio * (dpdx + mus * phi0) * (1e0 - x[i]) \
                * (hp - hc) ** 2 / (1e0 - hp)

    return uex, hc
# --------------------------------------------------------------
#
#            END OF THE FUNCTION
#
# --------------------------------------------------------------
