from __future__ import division, print_function, absolute_import
import numpy as np
from .actmodels_cy import nrtl_cy, rkter_nrtl_cy
from .actmodels_cy import dnrtl_cy, drkter_nrtl_cy


def nrtl_aux(X, tau, G):
    X = np.asarray(X, dtype=np.float64)
    lngama = nrtl_cy(X, tau, G)
    return lngama


def dnrtl_aux(X, tau, G):
    X = np.asarray(X, dtype=np.float64)
    lngama, dlngama = dnrtl_cy(X, tau, G)
    return lngama, dlngama


def nrtl(X, T, alpha, g, g1):
    '''
    NRTL activity coefficient model.

    Parameters
    ----------
    X: array like
        vector of molar fractions
    T: float
        absolute temperature in K
    g: array like
        matrix of energy interactions in K
    g1: array_like
        matrix of energy interactions in 1/K
    alpha: array_like
        aleatory factor.

    Notes
    -----
    tau = ((g + g1*T)/T)

    Returns
    -------
    lngama: array_like
        natural logarithm of activify coefficient
    '''
    X = np.asarray(X, dtype=np.float64)
    tau = g/T + g1
    G = np.exp(-alpha*tau)
    lngama = nrtl_cy(X, tau, G)

    return lngama


def nrtlter_aux(X, tau, G, D):
    X = np.asarray(X, dtype=np.float64)
    lngama = nrtl_cy(X, tau, G)

    xd = X*D
    lngama += rkter_nrtl_cy(X, xd)
    # lngama += nrtl(X, T, alpha, g, g1)
    return lngama


def dnrtlter_aux(X, tau, G, D):
    X = np.asarray(X, dtype=np.float64)

    xd = X*D
    lngamaD, dlngamaD = drkter_nrtl_cy(X, xd, D)
    lngama, dlngama = dnrtl_cy(X, tau, G)

    lngama += lngamaD
    dlngama += dlngamaD

    return lngama, dlngama


def nrtlter(X, T, alpha, g, g1, D):
    '''
    NRTL activity coefficient model.

    Parameters
    ----------
    X: array like
        vector of molar fractions
    T: float
        absolute temperature in K
    g: array like
        matrix of energy interactions in K
    g1: array_like
        matrix of energy interactions in 1/K
    alpha: array_like
        aleatory factor
    D : array_like
        ternary contribution parameters

    tau = ((g + g1*T)/T)

    Returns
    -------
    lngama: array_like
        natural logarithm of activify coefficient
    '''
    xd = X*D
    lngama = rkter_nrtl_cy(X, xd)
    lngama += nrtl(X, T, alpha, g, g1)
    return lngama


def dnrtl(X, T, alpha, g, g1):
    '''
    Derivatives of NRTL activity coefficient model.

    Parameters
    ----------
    X: array like
        vector of molar fractions
    T: float
        absolute temperature in K
    g: array like
        matrix of energy interactions in K
    g1: array_like
        matrix of energy interactions in 1/K
    alpha: array_like
        aleatory factor.

    Notes
    -----
    tau = ((g + g1*T)/T)

    Returns
    -------
    lngama: array_like
        natural logarithm of activify coefficient
    dlngama: array_like
        derivative of natural logarithm of activify coefficient
    '''
    X = np.asarray(X, dtype=np.float64)
    tau = g/T + g1
    G = np.exp(-alpha*tau)
    lngama, dlngama = dnrtl_cy(X, tau, G)

    return lngama, dlngama


def dnrtlter(X, T, alpha, g, g1, D):
    '''
    Derivatives of NRTL activity coefficient model with additional ternary
    contribution.

    Parameters
    ----------
    X: array like
        vector of molar fractions
    T: float
        absolute temperature in K
    g: array like
        matrix of energy interactions in K
    g1: array_like
        matrix of energy interactions in 1/K
    alpha: array_like
        aleatory factor
    D : array_like
        ternary contribution parameters

    tau = ((g + g1*T)/T)

    Returns
    -------
    lngama: array_like
        natural logarithm of activify coefficient
    dlngama: array_like
        derivative of natural logarithm of activify coefficient
    '''

    xd = X*D
    lngamaD, dlngamaD = drkter_nrtl_cy(X, xd, D)
    lngama, dlngama = dnrtl(X, T, alpha, g, g1)

    lngama += lngamaD
    dlngama += dlngamaD

    return lngama, dlngama
