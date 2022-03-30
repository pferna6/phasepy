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
Akl = 
Bkl = 
def DS(e):
    nums = e+1
    DSsum = 0
    for k in range(0,nums):
        for l in range(0,nums):
            DSsum += (al[0][k]-al[1][k])(al[0][l]-al[1][l])
    print(DSsum)
    
# alphaDS(2,2)
