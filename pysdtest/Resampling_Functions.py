""" 
Resampling Functions 

@Author: Kyungho Lee
Latest Update at July 19th 2021
"""

# Import Modules
import numpy as np

#%%
def subsampling(sample, subsize, nsub):
        
    """
        
    Subsampling
        
    Parameters
    ----------
    sample  : np.array (1-d)
    subsize : number of samples chosen by subsampling
    nsub    : number of subsample test statistics
        
    """
        
        
    # nsub : # of subsamples
    # subindex : subsample size x 1 x nsub
    subindex = (np.arange(0, subsize, 1)[:, None] + np.arange(0, nsub, 1)[:, None].T)[:, None]
    
    # subsample : subsample size x 1 x nsub x 1
    subsample = sample[subindex]
    return subsample

#%%
def bootstrap(sample, b, nbtsp):
            
    """
        
    Recentered Bootstrap
        
    Parameters
    ----------
    sample   : np.array (1-d)
    btspsize : number of samples chosen by bootstrapping
    nbtsp    : number of bootstrap test statistics
        
    """
        
    n = sample.shape[0]
    # nbtsp : # of bootsrap samples
    # btspindex : bootsrap size x 1 x nsamp
    btspindex = np.array([np.random.randint(n, size=b) for _ in np.arange(nbtsp)]).T[:,None,:]
    # btspsample : bootstrap size x 1 x nbtsp x 1 
    btspsample = sample[btspindex]
    return btspsample
        
#%%    
def paired_bootstrap(sample1, sample2, b, nbtsp):
        
    """
        
    Paired Bootstrap
        
    Parameters
    ----------
    sample1  : np.array (1-d)
    sample2  : np.array (1-d)
    btspsize : number of samples chosen by bootstrapping
    nbtsp    : number of bootstrap test statistics
        
    """
        
    n = sample1.shape[0]
    # nbtsp : # of bootsrap samples
    # btspindex : bootsrap size x 1 x nsamp
    btspindex = np.array([np.random.randint(n, size=b) for _ in np.arange(nbtsp)]).T[:,None,:]
    # btspsample : bootstrap size x 1 x nbtsp x 1 
    btspsample1 = sample1[btspindex]
    btspsample2 = sample2[btspindex]
        
    return btspsample1, btspsample2 
