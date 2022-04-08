# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:37:26 2022

@author: Pedro
"""
import numpy as np
from phasepy import mixture, component, preos
from phasepy.equilibrium import tpd_min, tpd_minimas

npropylbenzene = component(name='npropylbenzene', Tc=638.4, Pc=32, Zc=0.265, Vc=440, w=0.344,
                  GC={'ACH':5,'ACCH2':1,'CH2':1,'CH3':1})

undecane = component(name='undecane', Tc=639, Pc=19.8, Zc=.242, Vc=689, w=0.539,
                  GC={'CH2':9,'CH3':2})

nitrogen = component(name='nitrogen', Tc=126.2, Pc=33.98, Zc=.288, Vc=89, w=0.037,
                  GC={'N2':1})


mix = npropylbenzene+undecane+nitrogen
eos = preos(mix)

#combining individual components simply combines the values of the individual components in a file 

# mix.kij_cubic(Kij)
# eos = preos(mix, mixrule = 'qmr')
# T = 320.0
# P = 1.01
# z = np.array([0.5, 0.5])
# w = np.array([0.01, 0.99])
# print("Liquid molar fractions and TPD value:", tpd_min(w, z, T, P, eos, stateW='L', stateZ='L'))
# print("Vapor molar fractions and TPD value:", tpd_min(w, z, T, P, eos, stateW='V', stateZ='L'))

# x0 = np.array([0.23512692, 0.76487308])
# y0 = np.array([0.5, 0.5])
# flash(x0, y0, 'LV', Z, T, P, eos, full_output=True)

#%%
from Kij_modification_framework import kijCalc