# -*- coding: utf-8 -*-
"""
Created on Mon May 17 22:40:42 2021

Conversion of Antioine coefficients
10 to log
A = A/.434294
B = B/.434294
C = C


@author: Patrick
"""
"""UNITS:
    Temperature = K
    Pressure = bar
    volumes = cm**3/mol
    """

import numpy as np

from phasepy import component, mixture, preos
from phasepy.equilibrium import tpd_min, tpd_minimas, flash, bubblePy#, #dewPy
from PREOS_classes_PL import Jaubertkij
from collections import Counter


#%%
####### script #######
if __name__ == '__main__':
    # water = component(name='water', Tc=647.13, Pc=220.55, Zc=0.229, Vc=55.948,
    #                   Ant=[11.64785144, 3797.41566067, -46.77830444],
    #                       w=0.344861, GC={'H2O':1})
    # ethanol = component(name='ethanol', Tc=514.0, Pc=61.37, Zc=0.241, Vc=168.0,
    #                     Ant=[11.61809279, 3423.0259436, -56.48094263],
    #                         w=0.643558, GC={'CH3':1, 'CH2':1, 'OH(P)':1})
    npropylbenzene = component(name='n-Propylbenzene',Tc=638.4, Pc=32.0, Vc=440, Zc = .265, Mw= 120.2,
                               Ant=[9.384978,3433.073,-66.0],
                               w=.344 ,GC = {'ACH':5,'ACCH2':1,'CH2':1,'CH3':1}) #Unifac style, Jaubert does not have ACCH2
    undecane = component(name='Undecane',Tc = 639.0, Pc=19.8, Vc=689., Zc=.242, Mw=156.3,
                         Ant=[9.44439, 3620.766,-85.128],
                         w =.539, GC={'CH3':2, 'CH2':9}) #works for both Unifac and Jaubert
    nitrogen = component(name='Nitrogen',Tc = 126.2, Pc = 33.98, Vc = 89., Zc = 0.288, Mw = 28.01,
                         Ant = [8.603,609.382,-6.788],
                         w = 0.037, GC={'N2':1})  #Jaubert only.  unifac does not have N2.
    propane = component(name='Propane',Tc = 369.8, Pc = 42.455, Vc=200., Zc=0.277,Mw=44., #Vc and Zc come from engineering toolbox, not a proper reference
                    Ant=[10.44632, 2646.499, 24.906], #Come from NIST webook
                    w=0.152,GC={'CH3':2, 'CH2': 1})
    butane = component(name='Butane',Tc=425.2, Pc = 37.997, Vc=255., Zc = 0.274,Mw=58., #Vc and Zc come from engineering toolbox, not a proper reference
                       Ant = [10.02951,2706.875,-2.071],# NIST webbook Das, Reed, 1973
                       w = 0.193,GC={'CH3':2, 'CH2': 2})
    
    mix=mixture(npropylbenzene,undecane)
    # mix.add_component(propane)
    mix.add_component(nitrogen)
    # mix.add_component(butane)

#%%
    mix = Jaubertkij(mix,330)
    
    #Calculate binary interaction parameters from unifac.  We could redo this for Jaubert
    # mix.unifac()
    
    #also kij is formally temperature dependent in PREOS
    
    eos = preos(mix,'qmr') # Peng Robinson EOS using quadratic mixing rule
    

#%%

# if __name__ == '__main__':
#     import matplotlib.pyplot as plt
#     gases = [nitrogen,undecane,npropylbenzene]
#     names =[]
#     T_range = np.linspace(100,650,200)
#     vaporization_curve = np.ones((T_range.shape[0],3),)*np.nan
#     for i,gas in enumerate(gases):
#         names.append(gas.name)
#         for j in range(T_range.shape[0]):
#             if T_range[j] < gas.Tc:
#                 vaporization_curve[j,i]=gas.psat(T_range[j])
#         plt.scatter(gas.Tc,gas.Pc)
                
#     lines = plt.plot(T_range,vaporization_curve)
#     plt.legend(lines,names,fontsize = 15,bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
#            ncol=2, mode="expand", borderaxespad=0.)
#     plt.ylim([0,80])
#     plt.xlabel('Temperature / K',fontsize = 20)
#     plt.ylabel('Pressure / bar',fontsize = 20)
        
    
    
#%%
####### script #######
if __name__ == '__main__':

    T=623.05
    Pguess=72.15
    z=np.array([0.33333,0.33333,0.33333]) #global composition    
    
#     z=np.array([0.4,0.4,0.2]) #global composition
# #629.1
#     T=629.1
#     Pguess=48.1


    mix = Jaubertkij(mix,T)    
    eos = preos(mix,'qmr') # Peng Robinson EOS using quadratic mixing rule
    y0_guesses, tpd_values = tpd_minimas(10,z,T,Pguess,eos,'L','V') #guesses are self explanatory
    # jj = np.unique(np.array(y0_guesses),axis=0)
    x0_guesses, tpd_values2 = tpd_minimas(10,z,T,Pguess,eos,'V','L')
    
    #%%
    # ii = np.unique(np.array(x0_guesses),axis=0)
    
    # for i in range(ii.shape[0]):
    #     for j in range(jj.shape[0]):
    #         x0_guess = ii[i,:]
    #         y0_guess = jj[j,:]
    
    y0_guess = y0_guesses[np.argmax(tpd_values)]
    x0_guess = x0_guesses[np.argmin(tpd_values2)]
    #         # print(np.unique(x0_guesses),np.unique( tpd_values2))
    #         # x0_guess, tpd_value2 = tpd_min(W=w0,Z=z,T=T,P=Pguess,model=eos,stateW='L',stateZ='L')
    #         # print(x0_guess)
    x_result,y_result,fraction = flash(x_guess=x0_guess, y_guess = y0_guess, equilibrium='LV', Z=z,T=T, P=Pguess, model =eos)
    print(x_result-y_result) #This is the criteria for equilibra. For equilibra, the x and y component of a secific component must be approxietly equal. Cut off is .02? 
    print('The liquid phase mole fractions are', x_result)
    print('The vapor phase mole fractions are', y_result)
    print(fraction)
#%%

    # z_propane = np.linspace(0.2,.8,7)
    # x_result = np.zeros((7,2),)
    # y_result = np.zeros((7,2),)
    # P_result = np.zeros(7)
    # fraction = np.zeros(7)
    # for i in range(len(z_propane)):
    #     z= np.array([z_propane[i],1-z_propane[i]]) #global composition
    #     print(z)
    #     # for j in range(len

    #     # w0 = np.array([0.36, 0.63]) #initial value to try to find phase stability for , doesn't really matter, but you can try multiple 
    #     # y0_guess, tpd_value = tpd_min(W=w0,Z=z,T=T,P=Pguess, model=eos, stateW='V',stateZ='L')
    #     y0_guesses, tpd_values = tpd_minimas(5,z,T,Pguess,eos,'L','V')
    #     y0_guess = y0_guesses[np.argmax(tpd_values)]
    #     # print(y0_guess)
    #     x0_guesses, tpd_values = tpd_minimas(5,z,T,Pguess,eos,'V','L')
    #     x0_guess = x0_guesses[np.argmax(tpd_values)]
    #     # x0_guess, tpd_value2 = tpd_min(W=w0,Z=z,T=T,P=Pguess,model=eos,stateW='L',stateZ='L')
    #     # print(x0_guess)
    #     # x0_guesses, tpd_values = tpd_minimas(10, z, T, Pguess, eos, 'V', 'L')
    #     x_result[i,:],y_result[i,:],fraction[i] = flash(x_guess=x0_guess, y_guess = y0_guess, equilibrium='LV', Z=z,T=T, P=Pguess, model =eos)
    #     y_result[i,:], P_result[i] = bubblePy(y0_guess,Pguess,z,T,eos)

        