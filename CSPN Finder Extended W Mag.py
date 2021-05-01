# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 19:25:05 2020

@author: Isaac Radley

Attempting an analysis of a HASHPN X-match with GAIA EDR3 
to find best CSPN candidates 
"""

from PN_Source_analysis_withMag_RATIO import Analysis
import numpy as np
import astropy as ap
import pandas as pd


"""
Generating array from best match catalogue
"""
file=open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\EDR3 Files\GAIAArchiveEDR3_HashPN_Xmatch_ReducedCol_TEST.csv", 'r')
file_pd = open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\EDR3 Files\GAIAArchiveEDR3_HashPN_Xmatch_ReducedCol_TEST.csv", 'r')

"""
Retaining column names and finding column indices
"""
df = pd.read_csv(file_pd)
column = df.columns.tolist()


PnID_Index = df.columns.get_loc("idpnmain")
angDist_Index = df.columns.get_loc("dist")
GMag_Index = df.columns.get_loc('phot_g_mean_mag')

data = np.genfromtxt(file, delimiter=',', skip_header=1, dtype=float)

data[:,angDist_Index] *= 3600

ind1, ind2 = np.shape(data)

"""
Find the index at which the sources change as one HASH source will have many 
Gaia sources within the 60 arcsec range
"""
"""
Arrange data in order of idPNMain
"""
data = data[data[:,PnID_Index].argsort()]



Source_breaks = []

for i in range(0, len(data), 1):
    j = i+1
    
    if j == (len(data)):
        break
    
    elif data[i,PnID_Index] == data[j,PnID_Index]:
            
        continue
    else:
        Source_breaks.append(j)
             
Source_breaks = np.array(Source_breaks)        
    
"""
Building individual arrays from the source break index and the full data set
then saving each individual source as a file (probably not the most efficient)
"""

for i in range(0, len(Source_breaks)):
    q = Source_breaks[i]
    p = Source_breaks[i-1]
    
    if i == 0:
        
        y = data[0:q,:]
        Id = int(y[0,PnID_Index]) 
        y = pd.DataFrame(data=y, columns=column)
        y.to_csv(r'C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\PN_compTest_files\PN_ID_'+str(Id)+'.csv', index=False)
    else: 
        x = data[p:q,:]
        Id = int(x[0,PnID_Index])
        x = pd.DataFrame(data=x, columns=column)
        x.to_csv(r'C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\PN_compTest_files\PN_ID_'+str(Id)+'.csv', index=False)

"""
Analysing the individual Gaia sources around the HASHPN sources using 
my analysis fn.
"""
      
Final = []
column.append('Confidence')

for x in range(0, len(Source_breaks)):
    if x == 0 :
        sources = Analysis(int(data[0,PnID_Index]))
    else:
        Id = data[Source_breaks[x-1], PnID_Index]
        sources = Analysis(int(Id))
        
    if np.shape(sources) != (1, ind2 + 1):
        pass
    else:
        Final.append(sources)      

Final = np.array(Final)
Final = Final.reshape((len(Final),ind2 + 1))

Final = pd.DataFrame(data=Final, columns=column)
print(Final)

Final.to_csv(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\Best PN ratio files\BestPN_Presen_HighG.csv", index=False)