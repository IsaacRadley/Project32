# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 19:09:16 2021

@author: Isaac Radley
"""

import numpy as np
import pandas as pa
import math

def PNBackground(x):
    """ 
    Import GAIA EDR3 crossmatch with HASHPN file and retrieve location for 
    Parallax and Parallax/error
    """
    
    Gaia_file=open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\PN_compTest_files\PN_ID_" +str(x) +".csv", 'r')
    Gaia_file_pd = open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\PN_compTest_files\PN_ID_" +str(x) +".csv", 'r')
    Gaia_data = np.genfromtxt(Gaia_file, delimiter=',', skip_header=1)
    
    Gaia_df = pa.read_csv(Gaia_file_pd)
    column = Gaia_df.columns.tolist()
    
    Gaia_Parerr_index = Gaia_df.columns.get_loc("parallax_over_error")
    Gaia_Par_index = Gaia_df.columns.get_loc("parallax")
    Gaia_PnID_index = Gaia_df.columns.get_loc("idpnmain")
    
    
    """ 
    Import CSPN file and retrieve location for Parallax and Parallax/error
    """
    
    PN_file=open(r"INPUT_CSPN_CATALOGUE", 'r')
    PN_file_pd = open(r"INPUT_CSPN_CATALOGUE", 'r')
    PN_data = np.genfromtxt(PN_file, delimiter=',', skip_header=1)
    
    PN_df = pa.read_csv(PN_file_pd)
    
    PN_Par_index = PN_df.columns.get_loc("parallax")
    PN_Parerr_index = PN_df.columns.get_loc("parallax_over_error")
    PnID_index = PN_df.columns.get_loc("idpnmain")
    
    
    """
    Ensure Gaia array has correct shape 
    """
        
    a, b = np.shape(PN_data)
    
    if Gaia_data.ndim == 1:
        Gaia_data = np.reshape(Gaia_data, (-1,b-1))
    else:
        pass
    
    """
    Find Pn ID and whether or not that PN has a parallax value assinged to it,
    if not, pass
    """
    Pn_index = np.where(PN_data[:,PnID_index] == Gaia_data[0,Gaia_PnID_index])
    Pn_index = Pn_index[0]
    
    PN_background = []
    
    if np.shape(Pn_index) ==(0,):
        pass
    
    
    elif math.isnan(PN_data[Pn_index, PN_Par_index]) == True or math.isnan(PN_data[Pn_index, PN_Parerr_index]) == True :
        pass
    
    elif PN_data[Pn_index, PN_Parerr_index]<5 :
        pass
    
    else:
    
        """
        Remove all sources from the GAIA list with poor parallax measurements
        """
        
        good_parallax = []
        
        
        for i in range(0,len(Gaia_data)):
            if Gaia_data[i,Gaia_Parerr_index]< 5:
                pass
            else:
                y = Gaia_data[i,:]
                good_parallax.append(y)
                                          
                
        good_parallax = np.array(good_parallax)
        
        """
        By comparing parallaxes of the Best Pn catalogue and the GAIA good parallax
        catalogue, find sources at similar parallaxes.
        """
        
        
        for i in range(0, len(good_parallax)):
            if 0.8*PN_data[Pn_index,PN_Par_index] <= good_parallax[i,Gaia_Par_index] and good_parallax[i,Gaia_Par_index] <= 1.2*PN_data[Pn_index,PN_Par_index]:
                PN_background.append(good_parallax[i,:])        
            else: 
                pass
            
        PN_background = np.array(PN_background) 
        
    return PN_background

file=open(r"INPUT_CSPN_CATALOGUE", 'r')    
file_pd = open(r"INPUT_CSPN_CATALOGUE", 'r')

data = np.genfromtxt(file, delimiter=',', skip_header=1)

df = pa.read_csv(file_pd)
PnID_index = df.columns.get_loc("idpnmain")
column = df.columns.tolist()
column.remove('Confidence')

a,b = np.shape(data)

PN_Background = []
PN_Background = np.array(PN_Background)
PN_Background = np.reshape(PN_Background, (-1,b-1))

for x in range(len(data)):
    y = int(data[x,PnID_index])
    PN_Background_array = PNBackground(y)
    
    if np.shape(PN_Background_array) == (0,):
        pass
    else:
         PN_Background = np.concatenate((PN_Background,PN_Background_array))       

Final = pa.DataFrame(data=PN_Background, columns=column)
Final.to_csv(r"OUTPUT_BACKGROUND_STAR_CATALOGUE", index=False)