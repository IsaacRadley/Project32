# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:24:31 2020

@author: Isaac Radley

Actual analysis of PN source surroundings with ratio control of Brightness,
Colour and position
"""

import numpy as np
import pandas as pa
import math


def Analysis(x):
    file=open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\PN_compTest_files\PN_ID_" +str(x)+ ".csv", 'r')
    file_pd = open(r"C:\Users\kille\OneDrive\Desktop\Institute Of Astronomy\Project32 Characterising PNe\PN_compTest_files\PN_ID_" +str(x)+ ".csv", 'r')

    """
    Retaining column names and finding column indices
    """
    df = pa.read_csv(file_pd)
    column = df.columns.tolist()

    bp_rp_Index = df.columns.get_loc("bp_rp")
    angDist_Index = df.columns.get_loc("dist")
    PnRadius_Index = df.columns.get_loc("majdiam")
    GMag_Index = df.columns.get_loc('phot_g_mean_mag')
    
    data = np.genfromtxt(file, delimiter=',', skip_header=1)
    
    """
    Ensure the data has the right shape and dimensions
    """
    
    if data.ndim == (1):
        ind1 = np.shape(data)
        ind1 = ind1[0]
    else:
        ind1, ind2 = np.shape(data)

    Best_sources = []
    
    if np.shape(data) == (ind1,):
        data = np.reshape(data, (-1,ind1))
    else:
        pass
    """
    Finding the sources within half the PN Diam +2 arc secs by comparing to the
    angDist from the cross match
    """
    close_source = []
    
    for i in range(0, len(data), 1):
        data[i,PnRadius_Index] *= 0.5 
        """Half PN radius^^^"""
        
        if math.isnan(float(data[i,PnRadius_Index])) == True:
                
            if data[i,angDist_Index] < 10:
                    
                y = data[i,:]
                close_source.append(y)
            else:
                continue
        
        elif data[i,angDist_Index] < (2 + data[i,PnRadius_Index]):
            
            y = data[i,:]
            close_source.append(y)
        
        else:
            continue
                    
    close_source = np.array(close_source)
    
    
    """
    Analyse the astrometric and photometric data within the 
    pn + 2 cut above.
    
    """
    if np.shape(close_source) == (0,):
        pass
    else:
        close_nan_check = pa.DataFrame(data=close_source, columns=column)
        close_source[:,PnRadius_Index] *= 2
        
     
        """ If no colour or brightness then find closest to HASH """
        if close_nan_check['bp_rp'].isnull().all() == True and close_nan_check['phot_g_mean_mag'].isnull().all() == True:
            Min = np.nanmin(close_source[:,angDist_Index])
            index = np.where(close_source[:,angDist_Index] == Min )
            index = index[0]
            Best_sources = close_source[index,:]
        
        
        elif close_nan_check['bp_rp'].isnull().all() == True and close_nan_check['phot_g_mean_mag'].isnull().all() == False:
            Min = np.nanmin(close_source[:,GMag_Index])
            index = np.where(close_source[:,GMag_Index] == Min )
            index = index[0]
            Best_sources = close_source[index,:]
            """ If no colour value, find brightest source ^^^ """
        
        
        elif close_nan_check['bp_rp'].isnull().all() == False and close_nan_check['phot_g_mean_mag'].isnull().all() == True:
            Min = np.nanmin(close_source[:,bp_rp_Index])
            index = np.where(close_source[:,bp_rp_Index] == Min )
            index = index[0]
            Best_sources = close_source[index,:]
            """ If no brightness values, find bluest source ^^^ """
            
        else:
            
            """
            Arranges every star in close sources first by position then brightness then by colour and assigns
            a value based on the position within the list (highest number means lowest value
                                                           ie brighter or bluer)
            """
            length = list(range(0,len(close_source)))
            length.reverse()
            length = np.array(length)
            length = np.reshape(length, (-1,1))
            
            Ones = np.ones(np.shape(length))
            length = np.add(length,Ones)
        
            """Distance"""
            close_source = close_source[close_source[:,angDist_Index].argsort()]
            close_source = np.concatenate((close_source, length), axis=1)
            
            """Brightness"""
            close_source = close_source[close_source[:,GMag_Index].argsort()]
            close_source = np.concatenate((close_source, length), axis=1)
            
            """Colour"""
            close_source = close_source[close_source[:,bp_rp_Index].argsort()]
            close_source = np.append(close_source, length, axis=1)
        
            a, b = np.shape(close_source)
            
            """
            Ratio Weightings, B, G, P
            """
            close_source[:,b-1] *= 0.16
            close_source[:,b-2] *= 0.04
            close_source[:,b-3] *= 2.8
            
            Confidence = []

            for i in range(0, len(close_source), 1):
                y = close_source[i,(b-1)] + close_source[i,(b-2)] + close_source[i,(b-3)]
                Confidence.append(y)
            
            """
            Adds together each value to get a confidence estimate and finds the 
            highest confidence source 
            """
            Confidence = np.array(Confidence)
            Confidence_col = np.reshape(Confidence, (len(close_source), -1))
            
            close_source = np.append(close_source, Confidence_col, axis=1)
            
            close_source = np.delete(close_source, b-1, 1)
            close_source = np.delete(close_source, b-2, 1)
            close_source = np.delete(close_source, b-3, 1)
            
            a, b = np.shape(close_source)
            
            
            a = 3*a
            q = (100/a)
            
            close_source[:,b-1] *= q
            
            Max = np.nanmax(Confidence)
            index = np.where(Confidence == Max )
            index = np.array(index)
    
            if np.shape(index) == (1,1):
                index = index[0]
                Best_sources = close_source[index,:] 
        
            else:
                array = []
                array = np.array(array)
                array = np.reshape(array, (-1,b))
                a, b = np.shape(index)
                
                for i in range(0, b, 1):
                    x = index[0,i]
                    array = np.vstack((array, close_source[x,:]))
                
                Min = np.nanmin(array[:,bp_rp_Index])
                index = np.where(array[:,bp_rp_Index] == Min )
                index = index[0]
                Best_sources = array[index,:]
        
    return Best_sources
      


























