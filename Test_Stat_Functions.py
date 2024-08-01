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
        
    """
    
    Apply operator (as in Davidson and Duclos (2000)) for calculating (integrated) CDF

    
    Parameters
    ==============================
    sample: (1-dim) numpy array
        Observations          
    grid  : (ngrid x 1) numpy array
        Grid points
    s     : int
        Stochastic order (which determines the number of integration)

    return
    ==============================
        Results of applying operator which is a part of integration by parts in Davidson and Duclos (2000)

    """     
    
    if sample.ndim != 3:       
        sample = sample[:,None,None]
        
    # np.apply_along_axis(operator, axis, a, b)
    # apply operator (in here, op.lt) along axis
    # op.lt(a,b) <=> a < b
    # op.sub(a,b) <=> a - b

    return ((np.apply_along_axis(op.lt, 1, sample, grid) * (- np.apply_along_axis(op.sub, 1, sample, grid)) ** (
        s - 1)) / math.factorial(s - 1))

#%%    
def CDF(sample, grid, s):
    
    """

    Calculate (s-order integrated) CDF

    
    Parameters
    ==============================
    sample: (1-dim) numpy array
        Observations          
    grid  : (ngrid x 1) numpy array
        Grid points
    s     : int
        Stochastic order (which determines the number of integration)

    return
    ==============================
        CDF: (ngrid x 1) numpy array
            (integrated) CDF

    """       
    
    # We take average here.
    # if s = 1, D equals to ecdf
    CDF = operator(sample, grid, s).mean(0)
    
    return CDF

#%%
def set_grid(samples, ngrid):
    
    """
    
    Set grid points that are equally divided from the minimum to the maximum value of whole samples.
    
    Parameters
    ==============================
    samples : list
        A list contains samples, each sample is (1-d array) or Pandas Series.     
        
    ngrid   : int
        The number of grid points

    return
    ==============================
        (ngrid x 1) numpy array 
    
    """
    
    min_list = [np.min(samples[i]) for i in np.arange(len(samples))]
    max_list = [np.max(samples[i]) for i in np.arange(len(samples))]
    
    minx = np.min(np.min(min_list))
    maxx = np.max(np.max(max_list))     
    
    return np.linspace(minx,maxx, num=ngrid)