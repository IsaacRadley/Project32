# -*- coding: utf-8 -*-
"""

@author: Isaac Radley

Combining good parallax distances with 
Bailer Jones (2021) distance prior for poor parallax measurements
"""
import numpy as np
import pandas as pa

"""
Import Prior and PN files
"""
Prior_FileData = open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\Parallax_BailerJonesPrior_combiner\HASHxEDR3_DistPrior_Reduced.csv", 'r')
Prior_FileData_pd = open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\Parallax_BailerJonesPrior_combiner\HASHxEDR3_DistPrior_Reduced.csv", 'r')

PN_FileData = open(r"CSPN_INPUT_FILE", 'r')
PN_FileData_pd = open(r"CSPN_INPUT_FILE", 'r')

"""
Find column names and positions
"""
Prior_df = pa.read_csv(Prior_FileData_pd)
Prior_column = Prior_df.columns.tolist()

prior_sourceID = Prior_df.columns.get_loc("source_id")
prior_r_med_geo = Prior_df.columns.get_loc("r_med_geo")
prior_r_lo_geo = Prior_df.columns.get_loc("r_lo_geo")
prior_r_hi_geo = Prior_df.columns.get_loc("r_hi_geo")


PN_df = pa.read_csv(PN_FileData_pd)
PN_column = PN_df.columns.tolist()

PN_sourceID = PN_df.columns.get_loc("source_id")
PN_par_err = PN_df.columns.get_loc("parallax_over_error")
PN_parallax = PN_df.columns.get_loc("parallax")


"""
Create arrays from each file
"""

Prior_data = np.genfromtxt(Prior_FileData, delimiter=',', skip_header=1)
PN_data = np.genfromtxt(PN_FileData, delimiter=',', skip_header=1)

distances = []
prior_lo = []
prior_hi = []
Distance_flag = []

""" Count of how many surces with no data
"""
s = 0 

for x in range(0,len(PN_data)):
    """
    Sources with good parallax
    """
    if PN_data[x, PN_par_err] >= 5:
        dist = 1/PN_data[x, PN_parallax]
        distances.append(dist)
        prior_lo.append(0)
        prior_hi.append(0)
        Distance_flag.append(1)
        
    
    elif 1 <= PN_data[x, PN_par_err] < 5:
        """
        Sources requiring prior
        """
        source = PN_df.at[x,"source_id"]
        prior_index = Prior_df[Prior_df["source_id"]==source].index.values
        prior_index_int = prior_index[0]
        dist = Prior_data[prior_index_int, prior_r_med_geo]/1000
        distances.append(dist)
        prior_lo.append(Prior_data[prior_index_int, prior_r_lo_geo]/1000)
        prior_hi.append(Prior_data[prior_index_int, prior_r_hi_geo]/1000)
        Distance_flag.append(2)
        
    else:
        """
        Sources with no distance information
        """
        distances.append(0)
        prior_lo.append(0)
        prior_hi.append(0)
        s=s+1
        Distance_flag.append(0)
        
Distance_flag = np.array(Distance_flag)
Distance_flag = np.reshape(Distance_flag, (-1,1))

Distances = np.array(distances, dtype=float)
Distances = np.reshape(Distances, (-1,1))

Prior_lo = np.array(prior_lo, dtype=float)
Prior_lo = np.reshape(Prior_lo, (-1,1))

Prior_hi = np.array(prior_hi, dtype=float)
Prior_hi = np.reshape(Prior_hi, (-1,1))

PN_data = np.concatenate((PN_data,Distances), axis=1)
PN_data = np.concatenate((PN_data,Prior_lo), axis=1)
PN_data = np.concatenate((PN_data,Prior_hi), axis=1)
PN_data = np.concatenate((PN_data,Distance_flag), axis =1)

PN_column.append("Distance")
PN_column.append("r_lo_geo")
PN_column.append("r_hi_geo")
PN_column.append("Flag")

output = pa.DataFrame(data=PN_data, columns=PN_column)
print(output)

print(output.at[3,"source_id"])

output.to_csv(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\Parallax_BailerJonesPrior_combiner\BestPN75_DistanceCombination.csv", index=False)