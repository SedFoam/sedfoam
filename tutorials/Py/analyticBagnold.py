#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# --------------------------------------------------------------
#
#     AUTHOR : Julien CHAUCHAT                 DATE : 01/12/2019
#
# --------------------------------------------------------------
#
#     Function that solves the analytical solution for the
#       Bagnold profile
#
# --------------------------------------------------------------

"""
Created on Sun Dec  1 19:48:34 2019

@author: chauchat
"""


def analyticBagnold(nx, H, g, d, rho_p, rho_f, phi0, I0, Bphi, mus, mu2, beta):
    import numpy as np

    xex = np.linspace(0, H, nx)

    # Compute the analytical solution
    alphaex = np.ones(nx) * phi0

    I = I0 / ((mu2 - mus) / (np.tan(beta) - mus) - 1.0)
    alphaex = phi0 / (1 + Bphi * I) * np.ones(nx)
    uex = (
        2.0
        / 3.0
        * I
        * np.sqrt(alphaex * np.cos(beta))
        * (H**1.5 - (H - xex) ** (1.5))
    )

    #
    # dimensional form
    #
    U0 = np.sqrt(g * d)
    uex = uex * U0
    xex = xex * d

    pex = np.zeros(nx)
    pex = (rho_p - rho_f) * g * np.cos(beta) * alphaex * (H * d - xex)

    muIex = (mus + (mu2 - mus) / (I0 / I + 1)) * np.ones(nx)
    tauex = muIex * pex

    duexdz = np.zeros(nx)
    nuex = np.zeros(nx)
    for i in range(nx - 1):
        duexdz[i] = (uex[i + 1] - uex[i]) / (xex[i + 1] - xex[i])
        nuex[i] = (
            (mus + (mu2 - mus) / (I0 / I + 1))
            * pex[i]
            / (rho_p * (np.abs(duexdz[i]) + 1e-10))
        )

    return xex, I, alphaex, uex, pex, muIex, tauex, duexdz, nuex
