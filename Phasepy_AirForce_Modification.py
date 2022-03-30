# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:37:26 2022

@author: Pedro
"""

def doubleSummation(a,b):
    nums = a+1
    dums = b+1
    sum = 0
    for i in range(1,nums):
        for j in range(2,dums):
            sum +=i+j
    print(sum)
    
doubleSummation(3, 4)

#%%
def singleSummation(a):
    nums = a+1
    sum = 0
    for i in range(1,nums):
        sum += i
    print(sum)
    
singleSummation(3)

#%%Summation for beta 
def BetaSummation(a):
    nums = a+1
    bsum = 0
    z = [1,2,3,4]
    bi = [1,2,3,4]
    for i in range(1,nums):
        bsum += z[i]+bi[i]
    print(bsum)
    
BetaSummation(3)

#%%
# for i in range(1,3):
    # print(i)
Z = [1,2,3,4,4]
print(Z[0])   
#note, indexes for python starts at zero 

#%%Summation for alpha 
import numpy as np
import math as mt
Z = [1,2,3,4,5]
Y = [1,2,3,4,5]
k = np.array([[0,1,2],
              [1,0,3],
              [2,3,0]])

def alphaDS(a,b):
    nums = a+1
    dums = b+1
    alphasum = 0
    for i in range(0,nums):
        for j in range(0,dums):
            alphasum += Z[i]*Z[j]*mt.sqrt(Y[i]*Y[j])*(1-k[i][j])
    print(alphasum)
    
alphaDS(2,2)

#%% Summation for DS 
 
al = np.array([[2/3,1/3,0,0,0,0],
            [1/2,1/2,0,0,0,0]])
Akl = np.array([[0.0,74.81*10**6,261.5*10**6,396.7*10**6,32.94*10**6,8.579*10**6],
                [74.81*10**6,0.0,51.47*10**6,88.53*10**6,36.72*10**6,31.23*10**6],
                [261.5*10**6,51.47*10**6,0.0,-305.7*10**6,145.2*10**6,174.3*10**6],
                [396.7*10**6,88.53*10**6,-305.7*10**6,0.0,263.9*10**6,203.8*10**6],
                [32.94*10**6,36.72*10**6,145.2*10**6,263.9*10**6,0.0,13.04*10**6],
                [8.579*10**6,31.23*10**6,174.3*10**6,203.8*10**6,13.04*10**6,0.0]])

Bkl = np.array([[0.0,165.7*10**6,388.8*10**6,804.3*10**6,-35.00*10**6,-29.5*10**6],
                [165.7*10**6,0.0,79.61*10**6,315.00*10**6,108.4*10**6,84.76*10**6],
                [388.8*10**6,79.61*10**6,0.0,-250.8*10**6,301.6*10**6,352.1*10**6],
                [804.3*10**6,315.00*10**6,-250.8*10**6,0.0,531.5*10**6,203.8*10**6],
                [-35.00*10**6,108.4*10**6,301.6*10**6,531.5*10**6,0.0,6.863*10**6],
                [-29.5*10**6,84.76*10**6,352.1*10**6,203.8*10**6,6.863*10**6 ,0.0]])
T = 303.15
def DS(e):
    nums = e+1
    DSsum = 0
    for k in range(0,nums):
        for l in range(0,nums):
            DSsum += (al[0][k]-al[1][k])*(al[0][l]-al[1][l])
            # *Akl[0][1]*(298.15/T)**((Bkl[0][1]/Akl[0][1])-1)
    print(-1*DSsum)
    
DS(5)
