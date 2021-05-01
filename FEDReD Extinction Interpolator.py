# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 12:23:01 2021

@author: Isaac Radley
"""

import numpy as np
import pandas as pd

"""
Import FEDReD Extinction data
"""

file_ext=open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\FEDReD Extinction\A0_extinction.csv", 'r')
file_ext_pd = open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\FEDReD Extinction\A0_extinction.csv", 'r')

df_ext = pd.read_csv(file_ext_pd)

data_ext = np.genfromtxt(file_ext, delimiter=',', skip_header=1, dtype=float)

"""
Extinction upper and lower bounds
"""

file_A0Min=open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\FEDReD Extinction\A0_Min.csv", 'r')

data_A0Min = np.genfromtxt(file_A0Min, delimiter=',', skip_header=1, dtype=float)

file_A0Max=open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\FEDReD Extinction\A0_Max.csv", 'r')

data_A0Max = np.genfromtxt(file_A0Max, delimiter=',', skip_header=1, dtype=float)

"""
Import PN catalogue
"""

file_PN=open(r"INPUT_CSPN_CATALOGUE", 'r')
file_PN_pd = open(r"INPUT_CSPN_CATALOGUE", 'r')

df_PN = pd.read_csv(file_PN_pd)
column = df_PN.columns.tolist()

column.append('A0')

distance_Index = df_PN.columns.get_loc("Distance")
glon_Index = df_PN.columns.get_loc("l")

data_PN = np.genfromtxt(file_PN, delimiter=',', skip_header=1, dtype=float)



a,b = np.shape(data_PN)

"""
Find corresponding distance and Glon within each table
"""

longs = np.array(data_ext[0,:])
dists = np.array(data_ext[:,0])

A0 = []

non_values = np.empty(1,)
non_values[:] = np.nan

"""
Find the closest value in FEDReD to both Distance and galactic longitude of 
the desired CSPN.
"""

for i in range(0, len(data_PN), 1):
    glon_sep = []
    dist_sep = []
    
    for q in range(0, len(longs), 1):
        x = np.abs(data_PN[i,glon_Index] - longs[q])
        glon_sep.append(x)

    glon_sep = np.array(glon_sep)
    
    Min = np.nanmin(glon_sep)
    glon_arr_index = np.where(glon_sep == Min )
   
    for q in range(0, len(dists), 1):
        x = np.abs(data_PN[i,distance_Index] - dists[q])
        dist_sep.append(x)
       

    dist_sep = np.array(dist_sep)
    
    Min = np.nanmin(dist_sep)
    dist_arr_index = np.where(dist_sep == Min )
    
    y = data_ext[dist_arr_index[0],glon_arr_index[0]]
    
    """
    Dealing with no parallax value
    """
    
    if  np.shape(y) == (0,):
        y = non_values
   
        A0 = np.concatenate((A0,y))
        
    else:
        A0 = np.concatenate((A0,y))
        
A0 = np.reshape(A0, (-1,1))

data_PN = np.concatenate((data_PN, A0), axis=1)

Output = pd.DataFrame(data=data_PN, columns=column)

Output.to_csv(r"OUTPUT_CSPN_CATALOGUE", index=False)











