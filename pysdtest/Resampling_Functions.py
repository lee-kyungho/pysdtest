""" 
Resampling Functions 

@Author: Kyungho Lee
Latest Update at Aug 12th 2022
"""

# Import Modules
import numpy as np

#%%
def subsampling(sample, subsize, nsub):
        
    """
        
    Subsampling
        
    Parameters
    ==============================
    sample  : np.array (1-d)
    subsize : number of samples chosen by subsampling
    nsub    : number of subsample test statistics
        
    Returns
    ==============================
    subsample : (N x nsub) numpy array

    """
        
        
    # nsub : # of subsamples
    # subindex : subsample size x 1 x nsub
    subindex = (np.arange(0, subsize, 1)[:, None] + np.arange(0, nsub, 1)[:, None].T)[:, None]
    
    # subsample : subsample size x 1 x nsub x 1
    subsample = sample[subindex]
    return subsample

#%%
def bootstrap(sample, b, nboot):
            
    """
        
    Generate Bootstrap Sample
        
    Parameters
    ==============================
    sample   : np.array (1-d)
    btspsize : number of samples chosen by bootstrapping
    nboot    : number of bootstrap test statistics
        
    Returns
    ==============================
    btspsample : (N x nbtsp) numpy array
    
    """
        
    n = sample.shape[0]
    # nboot : # of bootsrap samples
    # btspindex : bootsrap size x 1 x nsamp
    btspindex = np.array([np.random.randint(n, size=b) for _ in np.arange(nboot)]).T[:,None,:]
    # btspsample : bootstrap size x 1 x nbtsp x 1 
    btspsample = sample[btspindex]
    return btspsample
        
#%%    
def paired_bootstrap(sample1, sample2, b, nboot):
        
    """
        
    Generate Paired Bootstrap Sample
        
    Parameters
    ==============================
    sample1  : (1-dim) numpy array
    sample2  : (1-dim) numpy array
    btspsize : number of samples chosen by bootstrapping
    nboot    : number of bootstrap test statistics

    ==============================
    btspsample : (N x nboot) numpy array
        
    """
        
    n = sample1.shape[0]
    # nboot : # of bootsrap samples
    # btspindex : bootsrap size x 1 x nsamp
    btspindex = np.array([np.random.randint(n, size=b) for _ in np.arange(nboot)]).T[:,None,:]
    # btspsample : bootstrap size x 1 x nboot x 1 
    btspsample1 = sample1[btspindex]
    btspsample2 = sample2[btspindex]
        
    return btspsample1, btspsample2 

