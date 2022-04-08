# -*- coding: utf-8 -*-
"""
Created on Sun May 23 18:53:00 2021
PREOS_classes_PL
Peng Robinson EOS classes (Molecule and Mixture)
Calculates EOS terms using Jaubert functional groups.

@author: Patrick
"""
import numpy as np
from pathlib import Path
from pandas import read_excel
from collections import Counter

#%%
R = 8.3145e-5  # universal gas constant, m3-bar/K-mol


#Lynch's path for the Jaubert database.
Jaubert_filename = Path(r"C:\Users\Pedro\phasepy\Jaubert_2013_PL1.xlsx") #r ebfore path means read
    

#%%
if Jaubert_filename.exists():  #Self explanatory
    A0 = read_excel(Jaubert_filename, sheet_name = 'interaction_A', index_col='Group', engine='openpyxl')*1e6 #load A factors. index_col is used as the row lobels of the Dataframe(two dimensional data structure (a table basically)). openpyxl supports newer excel formats
    A0.fillna(0, inplace=True) #replace 0 for unreadable. By setting true, the current dataframe is manipulated. IF fasle, a copy is created with the filled in zeros
    B0 = read_excel(Jaubert_filename, sheet_name = 'interaction_B', index_col='Group', engine='openpyxl')*1e6
    B0.fillna(0, inplace=True)
    A_jb = np.array(A0) #Turn into an array. Also replaces zeros in the upper triangles to have same reflected values as in lower triangle 
    B_jb = np.array(B0)
    
    for k in range(len(A0)):
        for l in range(len(A0)):
            if A_jb[k,l] == 0 and A_jb[l,k] !=0: #If value is 0 but the opposite diagonal is not 0
                A_jb[k,l] = A_jb[l,k] #take the A value that is not zero
                B_jb[k,l] = B_jb[l,k]
    
    #These are the mapping between unifac functional groups and Jaubert's groups.
    unifac_groups = read_excel(Jaubert_filename, 'unifac_translation',index_col='unifac',engine='openpyxl')
    unifac_groups.fillna(0, inplace=True)
    
    #These are the mapping of the Jaubert functional groups (backup.)
    jaubert_groups = read_excel(Jaubert_filename, 'Jaubert_translation',index_col='Jaubert_group',engine='openpyxl')
    jaubert_groups.fillna(0, inplace=True)
else:
    print('Jaubert database not found.  Cannot load kij')

#Lynch custom classes for PREOS solver.  
class Molecule:
    """
    Store molecule info here
    """
    def __init__(self, name, Tc, Pc, omega,n_g):
        """
        Pass parameters desribing molecules
        """
        #! name
        self.name = name
        #! Critical temperature (K)
        self.Tc = Tc
        #! Critical pressure (bar)
        self.Pc = Pc
        #! Accentric factor
        self.omega = omega
        #group coefficients
        self.n_g = n_g

    def a(self,T):
        """estimate PR a parameters for a component using critical point properties and accentricity"""
        Tr = T / self.Tc  # reduced temperature
        a0 = 0.457235 * (R*1e5)**2 * self.Tc**2 / (self.Pc*1e5)  #are pressure values transferred to pascal?
        if self.omega <=0.491:
            kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        else:
            kappa = 0.379642 + 1.48503 * self.omega - 0.164423 * self.omega**2 + 0.016666 * self.omega**3
        a = a0 * (1 + kappa * (1 - np.sqrt(Tr)))**2
        return(a)
    
    def b(self):
        """estimate PR b parameters for a component using critical point properties"""
        b = 0.0777961 * (R*1e5) * self.Tc / (self.Pc*1e5)  #Some sources have the constant as 0.7780
        return(b)
    
    def Psat(self,T):
        """estimate vapor pressure for component using critical point properties and accentricity"""
        return(np.exp((np.log(self.Pc) +  (np.log(10) * (7/3) * (1 + self.omega) * (1 - (self.Tc/T))))))
    
    def Ki(self,T,P):
        """Calculates Ki for a component using critical point properties and accentricity"""
        return(np.exp((np.log(self.Pc/P)) +  (np.log(10) * (7/3) * (1 + self.omega) * (1 - (self.Tc/T)))))
    
    def m(self):
        """Calculates the kappa value for a molecule from the accentricity"""
        if self.omega <=0.491:
            kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        else:
            kappa = 0.379642 + 1.48503 * self.omega - 0.164423 * self.omega**2 + 0.016666 * self.omega**3
        return(kappa)
    
    def Ng(self): #number of groups involved in calcultion 
        Ng = 0
        for i in self.n_g:
            Ng = Ng + i
        return(Ng)
    
    def print_params(self):
        """
        Print molecule parameters.
        """
        print("""Molecule: %s.
        \tCritical Temperature = %.1f K
        \tCritical Pressure = %.1f bar.
        \tAccentric factor = %f""" % (self.name, self.Tc, self.Pc, self.omega))
        
    def __repr__(self):
        string='\nclass Molecule containing:\n'
        for i,key in enumerate(self.__dict__): #enumerate takes collection or tuple and returns enumerate object. Assigns index to entries
            string+=str(key) +': '+ str(self.__dict__[key])+',\n' #dict is a dictionary or other mapping objective used to store. Combines entry from __init__ to general names.
        return(string)


class Mixture(object):
    """This is a container for molecules"""
    def __init__(self,mols,**kwargs): #**kwargs allows for dictionary to go through for key and value
        self.molecule_list = []
        if len(mols)!=0:
            self.molecule_list.extend(mols)
    
    def __getitem__(self,i):
        return self.molecule_list[i]
    
    def __repr__(self):
        string='\nclass Gas_Mixture containing:\n'
        for i,molecule in enumerate(self.molecule_list):
            string+=repr(self.molecule_list[i])+',\n'
        return(string)
    
    def add_molecule(self,molecule,**kwargs):
        "Add a molecule to the gas mixture"
        self.molecule_list.append(molecule)
        
    def alpha_ijkl(self,mol,ind):
        return(self.molecule_list[mol].n_g[ind]/self.molecule_list[mol].Ng())
    
    def a_mix(self,T,X):
        """calculate mixture a parameter"""
        a=0
        if type(X) == list:
            X=np.array(X)
        X= X/np.sum(X)
        for i in range(len(self.molecule_list)):
            for j in range(len(self.molecule_list)):
                aij = self.a_ij(T,i,j)
                a = a +aij*X[i]*X[j]
        return(a)
                
    def b_mix(self,X):
        """calculate mixture b parameter"""
        if type(X) == list:
            X=np.array(X)
        X= X/np.sum(X)
        b = 0
        for i in range(len(self.molecule_list)):
            b = b+ self.molecule_list[i].b()*X[i]
        return(b)
                
        
    def k_ij(self,T,ii=0,jj=1):
        """calculate binary interaction parameter using Jaubert group parameters"""
        ai = self.molecule_list[ii].a(T)
        bi = self.molecule_list[ii].b()
        aj = self.molecule_list[jj].a(T)
        bj = self.molecule_list[jj].b()
        
        k_1_sum = 0
        for k in range(len(self.molecule_list[0].n_g)):
            for l in range(len(self.molecule_list[0].n_g)):
                if A_jb[k,l]==0:
                    continue
                else:
                    k_1_sum = k_1_sum + (self.alpha_ijkl(ii,k)-self.alpha_ijkl(jj,k))*(self.alpha_ijkl(ii,l)-self.alpha_ijkl(jj,l))*A_jb[k,l]*(298.15/T)**(B_jb[k,l]/A_jb[k,l]-1)
        
        kij =  (-0.5*k_1_sum  -   (((np.sqrt(ai))/bi) - ((np.sqrt(aj))/bj))**2)/(2*((np.sqrt(ai*aj))/(bi*bj)))
        #Testing Prausnitz
        # kiijj = [[0,0.0033],[0.0033,0]]
        # kij = kiijj[ii][jj]
        
        return(kij)

    def k_ij_matrix(self,T=298):
        """Calculates the k_ij matrix for molecules using Jaubert group parameters.  Returns a symmetric matrix of length number of molecules"""
        number_of_molecules=len(self.molecule_list)
        kij_mat = np.zeros((number_of_molecules,number_of_molecules),)#Initialize matrix
        for i in range(number_of_molecules):
            for j in range(number_of_molecules):
                kij_mat[i,j]=self.k_ij(T,i,j) #call k_ij to solve for kij on that pair
        return(kij_mat)

    
    def a_ij(self,T,ii=0,jj=1):
        """calcualte a interaction coefficients"""
        aij = (self.molecule_list[ii].a(T) * self.molecule_list[jj].a(T))**0.5*(1-self.k_ij(T,ii,jj))
        return(aij)
        
    def fugacity_calc(self,T,P,z,x):
        """calculate fugacity for each component"""
        num_components = len(self.molecule_list)
        if type(x) == list:
            x=np.array(x)
        x= x/np.sum(x)
        a_mix = self.a_mix(T,x)
        A_mix = a_mix * P *1e5/ (R*1e5)**2 / T**2
        b_mix = self.b_mix(x)
        B_mix = b_mix * P *1e5/ (R*1e5) / T
        
        Lnfug = np.zeros(num_components)
        for i in range(num_components):
            averaged_Ai = 0
            for j in range(num_components):
                aij = self.a_ij(T,i,j)
                Aij = aij * P *1e5/ (R*1e5)**2 / T**2
                # print('A{0},{1}={2}'.format(i,j,Aij))
                averaged_Ai = averaged_Ai + Aij*x[j]
            # print(averaged_Ai)
            with np.errstate(all='ignore'):
                Lnfug[i] = -np.log(z - B_mix) + (z - 1.0) * self.molecule_list[i].b() / b_mix - A_mix / np.sqrt(8) / B_mix * (2.0*averaged_Ai/A_mix - self.molecule_list[i].b() / b_mix) *\
                    np.log((z + (1.0 + np.sqrt(2)) * B_mix) / (z + (1.0 - np.sqrt(2)) * B_mix))
        return(np.exp(Lnfug))


def Jaubertkij(phasepy_mix,T=298):
    """Calculates the Jaubert kij parameters from chemical functional groups loaded in phasepy format.
    Takes a phasepy mixture and a temperature (defaults 298K).
    Returns a phasepy mixture.  Slightly awkward but didn't want to modify the main phasepy mixture class
    """
    
    try: unifac_groups #Is unifac_groups in namespace?
    except:
        if Jaubert_filename.exists():
            
            unifac_groups = read_excel(Jaubert_filename, 'unifac_translation',index_col='unifac',engine='openpyxl')
            unifac_groups.fillna(0, inplace=True)
            jaubert_groups = read_excel(Jaubert_filename, 'Jaubert_translation',index_col='Jaubert_group',engine='openpyxl')
            jaubert_groups.fillna(0, inplace=True)
            A0 = read_excel(Jaubert_filename, 'interaction_A', index_col='Group', engine='openpyxl')*1e6
            A0.fillna(0, inplace=True)
            B0 = read_excel(Jaubert_filename, 'interaction_B', index_col='Group', engine='openpyxl')*1e6
            B0.fillna(0, inplace=True)
            
            A_jb = np.array(A0)
            B_jb = np.array(B0)
            
            for k in range(len(A0)):
                for l in range(len(A0)):
                    if A_jb[k,l] == 0 and A_jb[l,k] !=0:
                        A_jb[k,l] = A_jb[l,k]
                        B_jb[k,l] = B_jb[l,k]
        else:
            print('Jaubert database not found.  Cannot load kij')

    
    Ng = np.zeros((phasepy_mix.nc,len(A0)),int) #Initialize a 2d matrix of pure components x Jaubert functional groups.
#zeros entries: shape(integre or sequence of integer)
    puregc = phasepy_mix.GC # List of pure component functional groups.

    for i,gc in enumerate(puregc):
        subgroups = list(gc.keys()) # names
        counts = Counter() # initialize counter
        try:  #Map unifac functional groups to Jaubert groups
            jaubert_group = (unifac_groups.loc[subgroups, 'Jaubert'].values).astype(int)
        except: #map Jaubert groups.
            jaubert_group = (jaubert_groups.loc[subgroups, 'Jaubert_groupID'].values).astype(int)
        counts = Counter(dict(zip(jaubert_group, gc.values()))) #Count number of groups. Creates dictionary 
        try: #unifac groups may map to multiple Jaubert groups.
            jaubert_extra_group = (unifac_groups.loc[subgroups, 'extra'].values).astype(int)
            counts += Counter(dict(zip(jaubert_extra_group,gc.values())))    
            del(counts[0]) #Delete the 0s for when the extra groups do not add
        except:pass
        for j in range(len(A0)):
            Ng[i,j] = counts[j+1] #Populate the PL molecule functional group table.


    Jaubert_mix=Mixture([]) #Initialize an empty PL mixture
    for i in range(phasepy_mix.nc): # loop over phase_py mixture
        Jaubert_mix.add_molecule(Molecule(phasepy_mix.names[i],phasepy_mix.Tc[i],phasepy_mix.Pc[i],phasepy_mix.w[i],Ng[i,:])) #Put phasepy information into PL Molecule class
            
    kij = Jaubert_mix.k_ij_matrix(T) #Method in PL mixture class.
    phasepy_mix.kij_cubic(kij) #Put the kij's into the phasepy mixture.

    return phasepy_mix #Return updated phasepy_mixture.

