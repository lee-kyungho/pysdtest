"""
Functions for Computing Test Statistics

@Author: Kyungho Lee
Latest Update at July 19th 2021
"""

# Import Modules
import operator as op
import numpy as np
import math


#%%
def operator(sample, grid, s):
        
    # np.apply_along_axis(operator, axis, a, b)
    # apply operator (in here, op.lt) along axis
    # op.lt(a,b) <=> a < b
    # op.sub(a,b) <=> a - b
            
    return (np.apply_along_axis(op.lt, 1, sample, grid) * (- np.apply_along_axis(op.sub, 1, sample, grid)) ** (
        s - 1) / math.factorial(s - 1))

#%%    
def CDF(sample, grid, s):
        
    if sample.ndim != 3:       
        sample = sample[:,None,None]
        
    # We take average here.
    # if s = 1, D equals to ecdf
    CDF = operator(sample, grid, s).mean(0)
    
    return CDF

#%%
def set_grid(samples, ngrid):

    minx = np.min(np.min(samples))
    maxx = np.max(np.max(samples))     
    return np.linspace(minx,maxx, num=ngrid)