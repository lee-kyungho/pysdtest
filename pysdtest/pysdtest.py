"""

@author: Kyungho Lee (SNU)

Latest Update at July 5th 2021

Python Code for
Stochastic Dominance Test 

"""

# Import modules
import operator as op
import numpy as np
import matplotlib.pyplot as plt
import math
import random
import time


class test_sd :
    
    """

    Stochastic Dominance Test in Barett and Donald (2003 ECMA, hereafter BD), Linton, Maasoumi and Whang (2005 RES, hereafter LMW)

    Kolmogorov-Smirnov type Test Statistic via imposing the Least Favorable Case
    
    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0
    
    <=> 
    
    H0 : d_s <= 0 
    H1 : d_s > 0
    
    In here, we follow the notation in LMW

    Parameters
    ----------
    sample1 : np.array (1-d)
    sample2 : np.array (1-d)
                input samples should be 1-dim np.array
    
    grid    : int
                # of grid points
    s       : int
                Order of Stochastic Dominance
    b1      : int
                resampling size of the sample1.
    b2      : int
                resampling size of the sample2.
    resampling : str
                resampling should be one of 'subsampling' and 'bootstrap'  
                
                'subsampling' -> use sumbsampling method suggested by LMW
                'bootstrap'   -> use recentered bootstrap method as in BD and LMW
                
    nsamp   : int
                # of bootstrap statistics for bootstrap distribution
    
    """
    
    def __init__(self, sample1, sample2, ngrid, s, b1, b2, resampling, nsamp = 200) :
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.b1          = b1
        self.b2          = b2
        self.resampling  = resampling
        self.nsamp       = nsamp
            
        minx = min(sample1.min(0), sample2.min(0))
        maxx = max(sample1.max(0), sample2.max(0))     
        grid = np.linspace(minx,maxx, num=ngrid)
        
        self.grid = grid
        
    def testing(self) :
        
        """
        
        Returns
        -------
            None
            In here, we print the result
            
        test_stat   : float
                the value of the test statistic
        test_stat_b : numpy array
                values of the resampled test statistics
        pval        : float
                p-value of the test
        
        
        
        """
        
        sample1     = self.sample1
        sample2     = self.sample2
        ngrid       = self.ngrid
        s           = self.s
        b1          = self.b1
        b2          = self.b2
        resampling  = self.resampling
        nsamp       = self.nsamp
        
        
        start_time = time.time()
                
        # grid
        grid  = self.grid
        
        # Estimation
        test_stat   = self.T_N(sample1, sample2, grid, s)
        test_stat_b = self.resampled_stat(sample1, sample2, grid, s, b1, b2, resampling, nsamp)
        
        print(test_stat)
        print(test_stat_b)
        
        pval = (test_stat_b >= test_stat).mean(0) 
        
        if s == 1:
            Torder = 'First Order SD'
        elif s == 2:
            Torder = 'Second Order SD'
        elif s == 3:
            Torder = 'Third Order SD'
        else:
            Torder = str(s) + 'th order SD'
    
    

        print('\n#-------------------------------------------#')    
        print('Testing Stochastic Dominance',
          '\n\n* H0 : sample1 ', Torder, 'sample2')
        
        print('* Resampling method \t:', resampling)
        print('\n')
        
        print('* # (sample1) \t\t = %6d' % sample1.shape[0],
          '\n* # (sample2) \t\t = %6d\n' % sample2.shape[0])
        print('* # ('+ resampling + '1) \t = %6d' % b1,
          '\n* # ('+ resampling + '2) \t = %6d\n' % b2)
        print('#-------------------------------------------#\n')    
        print('* SD order \t\t = %6d' % s,
          '\n* # of grid points \t = %6d\n' % ngrid)
        print('#-------------------------------------------#\n')    
        print('* Test Result *\n')    

        print('* Test statistic \t = %5.4f' % test_stat)
        print('* p-value \t\t = %5.4f\n' % pval)
        print('#-------------------------------------------#')    
        et = time.time() - start_time
        print('\n* Time elapsed : %5.2f Sec' % et)            
    
        self.result = {'test_stat'   : test_stat[0,0] ,
                   'test_stat_b' : np.squeeze(test_stat_b),
                   'pval'        : pval[0]
                   }
        
    def operator(self, sample, grid, s):
        
        # np.apply_along_axis(operator, axis, a, b)
        # apply operator (in here, op.lt) along axis
        # op.lt(a,b) <=> a < b
        # op.sub(a,b) <=> a - b
            
        return (np.apply_along_axis(op.lt, 1, sample, grid) * (- np.apply_along_axis(op.sub, 1, sample, grid)) ** (
            s - 1) / math.factorial(s - 1))
    
    def D(self, sample, grid, s):
        
        if sample.ndim != 3:       
            sample = sample[:,None,None]
        
        # We take average here.
        # if s = 1, D equals to ecdf
        
        operator = self.operator
        D = operator(sample, grid, s).mean(0)
        
        return D 
    
    def T_N(self, sample1, sample2, grid, s):
            
        """
            
        test statsitics
        
        d_s in LMW
        S_jf in BD
            
        """
        
        D = self.D
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = D(sample1, grid, s) - D(sample2, grid, s)
        
        if D_s.ndim == 1:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)
        
        else:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)[:,None]
        
    def subsampling(self, sample, subsize, nsub):
        
        """
        
        Subsampling
        
        Parameters
        ----------
        sample  : np.array (1-d)
        subsize : number of samples chosen by subsampling
        nsum    : number of subsample test statistics
        
        """
        
        
        # nsub : # of subsamples
        # subindex : subsample size x 1 x nsub
        subindex = (np.arange(0, subsize, 1)[:, None] + np.arange(0, nsub, 1)[:, None].T)[:, None]
    
        # subsample : subsample size x 1 x nsub x 1
        subsample = sample[subindex]
        return subsample
    
    def bootstrap(self, sample, btspsize, nbtsp):
        
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample = sample[btspindex]
        return btspsample
        
    
    def paired_bootstrap(self, sample1, sample2, btspsize, nbtsp):
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample1 = sample1[btspindex]
        btspsample2 = sample2[btspindex]
        
        return btspsample1, btspsample2 
    
    def resampled_stat(self, sample1, sample2, grid, s, b1, b2, resampling, nsamp):
        
        D = self.D
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
            
        if resampling == 'subsampling' :
            
            subsampling = self.subsampling

            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            lbd = n2 / (n1 + n2)
            b = lbd * b1
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = D(subsample1, grid, s) - D(subsample2, grid, s) # ngrid x N-b+1
            
            # take supremum over support
            return (b) ** (0.5) * (D_sub).max(0)[:, None]
                
        elif resampling == 'bootstrap' :
            bootstrap = self.bootstrap

            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nsamp)
            btspsample2 = bootstrap(sample2, b2, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (D_btsp).max(0)[:, None]
        
        elif resampling == "paired_bootstrap" :
            
            paired_bootstrap = self.paired_bootstrap
            
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (D_btsp).max(0)[:, None]
    
        else :
            print("check resampling parameter")
      
        
    def plot_D(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        
        Parameters
        ==========
        save  : boolean
            if True, save figure
        titel : str
            name of the figure to save
        label 1,2 : str
            the name of legend
        xlabel : str
            the label of x-axis.
        
        
        """
        
        
        sample1 = self.sample1
        sample2 = self.sample2
        
        grid    = self.grid
        s       = self.s
        D       = self.D

        D1 = D(sample1,grid,s)
        D2 = D(sample2,grid,s)
        
        # if s == 1 :
        #     ylabel = "Cumulative distribution function"
        # elif s > 1 :
        #     ylabel = str(s) + "-time cumulated distribution function"
            
        plt.figure()
        plt.plot(grid, D1, label = label1)
        plt.plot(grid, D2, label = label2, ls = 'dashed')
        plt.xlabel(xlabel)
        # plt.ylabel(ylabel)
        plt.legend()
        
        if save == True:
            plt.savefig(title)

        plt.show()
        
        
    
        
class test_sd_contact :
  
    """

    Stochastic Dominance Testing in Linton, Song and Whang (2010 J.Econometrics, LSW)

    L-2 type Test Statistic 
    Contact-set estimation is used for constructing resampled statistics
    
    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0
    
    <=> 
    
    H0 : d_s <= 0 
    H1 : d_s > 0
    
    In here, we follow the notation in LSW

    Parameters
    ----------
    sample1 : np.array (1-d)
    sample2 : np.array (1-d)
                input samples should be 1-dim np.array
    
    grid    : int
                # of grid points
    s       : int
                Order of Stochastic Dominance
    b1      : int
                resampling size of the sample1.
    b2      : int
                resampling size of the sample2.
    resampling : str
                resampling should be one of 'subsampling', 'bootstrap' and 'paired_bootstrap'
                
                'subsampling' -> use sumbsampling method suggested by LMW
                'bootstrap'   -> use recentered bootstrap method as in BD and LMW
                
    nsamp   : int
                # of bootstraping (default value: 200)
    
    c       : float
                Tuning parameter for contact-set estimation (default value: 0.75)
    
    """
    
    def __init__(self, sample1, sample2, ngrid, s, b1, b2, resampling, nsamp = 200, c = 0.75) :
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.b1          = b1
        self.b2          = b2
        self.nsamp       = nsamp
        self.n1          = sample1.shape[0]   
        self.n2          = sample2.shape[0]   
        self.resampling  = resampling


        # setting grid
        minx = min(sample1.min(0), sample2.min(0))
        maxx = max(sample1.max(0), sample2.max(0))     
        grid = np.linspace(minx,maxx, num=ngrid)
        
        self.grid = grid
        
        
        # Tuning parameter
        self.c           = c
        N = (self.n1 + self.n2) / 2
        self.N           = N
        # sequence c_n for the contact set estimation
        self.seq_c       = c*np.log(np.log(N))/np.sqrt(N)
    
    def testing(self) :
        
        """
        
        Returns
        -------
            None
            In here, we print the result
            
        test_stat   : float
                the value of the test statistic
        test_stat_b : numpy array
                values of the resampled test statistics
        pval        : float
                p-value of the test
        
        
        
        """
        
        sample1     = self.sample1
        sample2     = self.sample2
        ngrid       = self.ngrid
        s           = self.s
        b1          = self.b1
        b2          = self.b2
        resampling  = self.resampling
        nsamp       = self.nsamp
        c           = self.c
        
        start_time = time.time()
                
        # grid
        grid  = self.grid
        
        # Estimation
        test_stat   = self.T_N(sample1, sample2, grid, s)
        test_stat_b = self.resampled_stat(sample1, sample2, grid, s, b1, b2, resampling, nsamp)
        
        pval = (test_stat_b >= test_stat).mean(0) 
        
        if s == 1:
            Torder = 'First Order SDes'
        elif s == 2:
            Torder = 'Second Order SDes'
        elif s == 3:
            Torder = 'Third Order SDes'
        else:
            Torder = str(s) + 'th order SDes'
    
    
        print('\n#-------------------------------------------#')    
        print('Testing Stochastic Dominance by the Contact Set Estimation',
          '\n\n* H0 : sample1 ', Torder, 'sample2')
        print('Linton, Song, and Whang (2010, LSW)')    
        
        print('* Resampling method \t:', resampling)
        print('\n')
        
        print('* # (sample1) \t\t = %6d' % sample1.shape[0],
          '\n* # (sample2) \t\t = %6d\n' % sample2.shape[0])
        print('* # ('+ resampling + '1) \t = %6d' % b1,
          '\n* # ('+ resampling + '2) \t = %6d\n' % b2)
        print('#-------------------------------------------#\n')    
        print('* SD order \t\t = %6d' % s,
          '\n* # of grid points \t = %6d\n' % ngrid)
        print("Tuning parameter c = %6d" % c)
        print('#-------------------------------------------#\n')    
        print('* Test Result *\n')    

        print('* Test statistic \t = %5.4f' % test_stat)
        print('* p-value \t\t = %5.4f\n' % pval)
        print('#-------------------------------------------#')    
        et = time.time() - start_time
        print('\n* Time elapsed : %5.2f Sec' % et)            
    
        self.result = {'test_stat'   : test_stat[0,0] ,
                   'test_stat_b' : np.squeeze(test_stat_b),
                   'pval'        : pval[0]
                   }
        
        
    def operator(self, sample, grid, s):
        
        # np.apply_along_axis(operator, axis, a, b)
        # apply operator (in here, op.lt) along axis
        # op.lt(a,b) <=> a < b
        # op.sub(a,b) <=> a - b
            
        return (np.apply_along_axis(op.lt, 1, sample, grid) * (- np.apply_along_axis(op.sub, 1, sample, grid)) ** (
            s - 1) / math.factorial(s - 1))
    
    def D(self, sample, grid, s):
        
        if sample.ndim != 3:       
            sample = sample[:,None,None]
        
        # We take average here.
        # if s = 1, D equals to ecdf
        
        operator = self.operator
        D = operator(sample, grid, s).mean(0)
        
        return D 
    
    def T_N(self, sample1, sample2, grid, s):
            
        """
            
        test statsitics
        
        d_s in LMW
        S_jf in BD
            
        """
        
        D = self.D
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = D(sample1, grid, s) - D(sample2, grid, s)
        
        self.D_s  = D_s
        
        
        if D_s.ndim == 1:
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_s>0) * D_s)**(2),axis=0)
        
        else:
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_s>0) * D_s)**(2),axis=0)[:,None]
        
    def subsampling(self, sample, subsize, nsub):
        
        """
        
        Subsampling
        
        Parameters
        ----------
        sample  : np.array (1-d)
        subsize : number of samples chosen by subsampling
        nsum    : number of subsample test statistics
        
        """
        
        
        # nsub : # of subsamples
        # subindex : subsample size x 1 x nsub
        subindex = (np.arange(0, subsize, 1)[:, None] + np.arange(0, nsub, 1)[:, None].T)[:, None]
    
        # subsample : subsample size x 1 x nsub x 1
        subsample = sample[subindex]
        return subsample
    
    def bootstrap(self, sample, btspsize, nbtsp):
        
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample = sample[btspindex]
        return btspsample
        
    
    def paired_bootstrap(self, sample1, sample2, btspsize, nbtsp):
        
        """
        
        Paired Bootstrap
        
        Bootstrapping (X_1, X_2) together to allow dependence.
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample1 = sample1[btspindex]
        btspsample2 = sample2[btspindex]
        
        return btspsample1, btspsample2 
    
    def contact_set_estimation(self):
    
        """
        
        Contact-set estimation by Linton, Song and Whang (2010)
        This estimation is condcuted by using original samples.
        
        Parameters
        ----------
            None
            
        
        Results
        ----------
        the following variables are generated 
        
        contact_set : numpy array (ngrid x 1)
        ct_num      : int 
        
        """
        
        seq_c = self.seq_c
        D_s   = self.D_s
        
        # estimation
        # pick indices for the 'contact set'
        contact_set = np.abs(D_s) < seq_c
        
        # to identify whether the contact set is empty
        ct_num      = np.sum(contact_set)
        
        self.ct_num = ct_num
        
        if ct_num == 0:
            self.contact_set = np.ones(D_s.shape[0])
        
        self.contact_set = contact_set
        
    
    def resampled_stat(self, sample1, sample2, grid, s, b1, b2, resampling, nsamp):
        
        D                      = self.D
        contact_set_estimation = self.contact_set_estimation
        n1                     = self.n1
        n2                     = self.n2
        ngrid                  = self.ngrid 
        # contact-set estimation using the whole sample
        contact_set_estimation()        
        contact_set = self.contact_set
        ct_num      = self.ct_num
        
        
        if resampling == 'subsampling' :
            
            subsampling = self.subsampling

            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            lbd = n2 / (n1 + n2)
            b = lbd * b1
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = D(subsample1, grid, s) - D(subsample2, grid, s) # ngrid x N-b+1
            
            # generate contact_set_rep for bootstrap sample
            # contact_set_rep is ngrid x nsamp dimension
            
            contact_set_rep = np.repeat(contact_set,nsub) # (ngrid * nsub) x 1
            contact_set_rep = np.reshape(contact_set_rep, [ngrid, nsub]) # ngrid x nsub

            D_sub_contact = D_sub[contact_set_rep] # (ct_num * nsub) x 1
            D_sub_contact = np.reshape(D_sub_contact,[ct_num,nsub]) # ct_num x nsub
            
            
            # take integral on contact-set
            return (b) * np.trapz(((D_sub_contact>0) * D_sub_contact)**(2),axis=0)[:, None]
                
        elif resampling == 'bootstrap' :
            bootstrap = self.bootstrap

            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nsamp)
            btspsample2 = bootstrap(sample2, b2, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
        
            # generate contact_set_rep for bootstrap sample
            # contact_set_rep is ngrid x nsamp dimension
            
            contact_set_rep = np.repeat(contact_set,nsamp) # (ngrid * nsamp) x 1
            contact_set_rep = np.reshape(contact_set_rep, [ngrid, nsamp]) # ngrid x nsamp

            D_btsp_contact = D_btsp[contact_set_rep]
            D_btsp_contact = np.reshape(D_btsp_contact,[ct_num,nsamp])
            
            # take supremum over support
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_btsp_contact>0) * D_btsp_contact)**(2),axis=0)[:,None]
        
        elif resampling == "paired_bootstrap" :
            
            paired_bootstrap = self.paired_bootstrap
            
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # generate contact_set_rep for bootstrap sample
            # contact_set_rep is ngrid x nsamp dimension
            
            contact_set_rep = np.repeat(contact_set,nsamp) # (ngrid * nsamp) x 1
            contact_set_rep = np.reshape(contact_set_rep, [ngrid, nsamp]) # ngrid x nsamp

            D_btsp_contact = D_btsp[contact_set_rep]
            D_btsp_contact = np.reshape(D_btsp_contact,[ct_num,nsamp])
                        
            # take supremum over support
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_btsp_contact>0) * D_btsp_contact)**(2),axis=0)[:,None]
    
        else :
            print("check resampling parameter")
      
        
    def plot_D(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        
        Parameters
        ==========
        save  : boolean
            if True, save figure
        titel : str
            name of the figure to save
        label 1,2 : str
            the name of legend
        xlabel : str
            the label of x-axis.
        
        
        """
        
        
        sample1 = self.sample1
        sample2 = self.sample2
        
        grid    = self.grid
        s       = self.s
        D       = self.D

        D1 = D(sample1,grid,s)
        D2 = D(sample2,grid,s)
        
        # if s == 1 :
        #     ylabel = "Cumulative distribution function"
        # elif s > 1 :
        #     ylabel = str(s) + "-time cumulated distribution function"
            
        plt.figure()
        plt.plot(grid, D1, label = label1)
        plt.plot(grid, D2, label = label2, ls = 'dashed')
        plt.xlabel(xlabel)
        # plt.ylabel(ylabel)
        plt.legend()
        
        if save == True:
            plt.savefig(title)

        plt.show()
        
        
################# Donald and Hsu ###########################
        
class test_sd_SR :
  
    """

    Stochastic Dominance Testing in Donald and Hsu (2014, Econemtrics Review)

    KS-type Test Statistic 
    Selective recentering approach
    
    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0
    
    <=> 
    
    H0 : d_s <= 0 
    H1 : d_s > 0
    
    Parameters
    ----------
    sample1 : np.array (1-d)
    sample2 : np.array (1-d)
                input samples should be 1-dim np.array
    
    grid    : int
                # of grid points
    s       : int
                Order of Stochastic Dominance
    b1      : int
                resampling size of the sample1.
    b2      : int
                resampling size of the sample2.
    resampling : str
                resampling should be one of 'subsampling', 'bootstrap', and 'paired_bootstrap'   
                
                'subsampling' -> use sumbsampling method as in LMW
                'bootstrap'   -> use recentered bootstrap method as in BD and LMW
                'paired_bootstrap' -> use recentered bootstrap method by resampling a pair (X,Y) to allow dependency of paired observations.
                
    nsamp   : int
                # of bootstrap statistics for bootstrap distribution   
    a       : float
                Tuning parameter for selective-recentering approach. (default value: 0.1)
    eta     : float
                Tuning parameter for selective-recentering approach. (default value: 10^(-6))

    """
    
    def __init__(self, sample1, sample2, ngrid, s, b1, b2, resampling, nsamp = 200, a = 0.1, eta = 10**(-6)):
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.b1          = b1
        self.b2          = b2
        self.nsamp       = nsamp
        
        n1 = sample1.shape[0] 
        n2 = sample2.shape[0] 
        
        self.n1          = sample1.shape[0]   
        self.n2          = sample2.shape[0]   
        self.resampling  = resampling
        self.a           = a

        # setting grid
        minx = min(sample1.min(0), sample2.min(0))
        maxx = max(sample1.max(0), sample2.max(0))     
        grid = np.linspace(minx,maxx, num=ngrid)
        
        self.grid = grid
        
        seq_a = -a* np.sqrt(np.log (np.log(n1+n2)))
        # Tuning parameter
        self.seq_a           = seq_a
        self.eta             = eta
        
    def testing(self) :
        
        """
        
        Returns
        -------
            None
            In here, we print the result
            
        test_stat   : float
                the value of the test statistic
        test_stat_b : numpy array
                values of the resampled test statistics
        pval        : float
                p-value of the test
        
        
        
        """
        
        sample1     = self.sample1
        sample2     = self.sample2
        ngrid       = self.ngrid
        s           = self.s
        b1          = self.b1
        b2          = self.b2
        resampling  = self.resampling
        nsamp       = self.nsamp
        a           = self.a
        eta         = self.eta
        
        
        start_time = time.time()
                
        # grid
        grid  = self.grid
        
        # Estimation
        test_stat   = self.T_N(sample1, sample2, grid, s)
        self.test_stat = test_stat
        
        test_stat_b = self.resampled_stat(sample1, sample2, grid, s, b1, b2, resampling, nsamp)
        test_stat_b[test_stat_b <= 10**(-6)] = 10^(-6)
        
        pval = (test_stat_b >= test_stat).mean(0) 
        
        if s == 1:
            Torder = 'First Order SD'
        elif s == 2:
            Torder = 'Second Order SD'
        elif s == 3:
            Torder = 'Third Order SD'
        else:
            Torder = str(s) + 'th order SD'
    
    

        print('\n#-------------------------------------------#')    
        print('Testing Stochastic Dominance by Selective Recentering Approach (Donald and Hsu 2014)',
          '\n\n* H0 : sample1 ', Torder, 'sample2')
        
        print('* Resampling method \t:', resampling)
        print('\n')
        
        print('* # (sample1) \t\t = %6d' % sample1.shape[0],
          '\n* # (sample2) \t\t = %6d\n' % sample2.shape[0])
        print('* # ('+ resampling + '1) \t = %6d' % b1,
          '\n* # ('+ resampling + '2) \t = %6d\n' % b2)
        print('#-------------------------------------------#\n')    
        print('* SD order \t\t = %6d' % s,
          '\n* # of grid points \t = %6d\n' % ngrid)
        print('# Tuning paremeters -------------')
        print('# a = %6d' % a)
        print('# eta = %6d' % eta)
        print('#-------------------------------------------#\n')    
        print('* Test Result *\n')    

        print('* Test statistic \t = %5.4f' % test_stat)
        print('* p-value \t\t = %5.4f\n' % pval)
        print('#-------------------------------------------#')    
        et = time.time() - start_time
        print('\n* Time elapsed : %5.2f Sec' % et)            
    
        self.result = {'test_stat'   : test_stat[0,0] ,
                   'test_stat_b' : np.squeeze(test_stat_b),
                   'pval'        : pval[0]
                   }
        
    def operator(self, sample, grid, s):
        
        # np.apply_along_axis(operator, axis, a, b)
        # apply operator (in here, op.lt) along axis
        # op.lt(a,b) <=> a < b
        # op.sub(a,b) <=> a - b
            
        return (np.apply_along_axis(op.lt, 1, sample, grid) * (- np.apply_along_axis(op.sub, 1, sample, grid)) ** (
            s - 1) / math.factorial(s - 1))
    
    def D(self, sample, grid, s):
        
        if sample.ndim != 3:       
            sample = sample[:,None,None]
        
        # We take average here.
        # if s = 1, D equals to ecdf
        
        operator = self.operator
        D = operator(sample, grid, s).mean(0)
        
        return D 
    
    def T_N(self, sample1, sample2, grid, s):
            
        """
            
        test statsitics
        
        S_jf in BD
            
        """
        
        D = self.D
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = D(sample1, grid, s) - D(sample2, grid, s)
        self.D_s = D_s
        
        
        if D_s.ndim == 1:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)
        
        else:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)[:,None]

        
    def subsampling(self, sample, subsize, nsub):
        
        """
        
        Subsampling
        
        Parameters
        ----------
        sample  : np.array (1-d)
        subsize : number of samples chosen by subsampling
        nsum    : number of subsample test statistics
        
        """
        
        
        # nsub : # of subsamples
        # subindex : subsample size x 1 x nsub
        subindex = (np.arange(0, subsize, 1)[:, None] + np.arange(0, nsub, 1)[:, None].T)[:, None]
    
        # subsample : subsample size x 1 x nsub x 1
        subsample = sample[subindex]
        return subsample
    
    def bootstrap(self, sample, btspsize, nbtsp):
        
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample = sample[btspindex]
        return btspsample
        
    
    def paired_bootstrap(self, sample1, sample2, btspsize, nbtsp):
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample1 = sample1[btspindex]
        btspsample2 = sample2[btspindex]
        
        return btspsample1, btspsample2 
    
    def selective_recentering(self):
    
        """
        
        Selective Recentering by Donald and Hsu (2014)
        This estimation is condcuted by using original samples.
        
        Parameters
        ----------
            None
            
        
        Results
        ----------
        the following variables are generated 
        
        selected_set : numpy array (ngrid x 1)
        sr_num      : int 
        
        """
        
        seq_a = self.seq_a
        D_s   = self.D_s
        n1    = self.n1
        n2    = self.n2
        
        # estimation
        # pick indices for the 'contact set'
        selected_set = D_s < seq_a
        
        # to identify whether the contact set is empty
        sl_num      = np.sum(selected_set)
        
        self.sl_num = sl_num        
        self.selected_set = selected_set
        
        return (n1 * n2 /(n1 + n2)) ** (0.5) * (D_s) * selected_set
        
    
    
    def resampled_stat(self, sample1, sample2, grid, s, b1, b2, resampling, nsamp):
        
        D = self.D
        selective_recentering = self.selective_recentering
        
        recentering_function = selective_recentering()
        
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
            
        if resampling == 'subsampling' :
            
            subsampling = self.subsampling

            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            lbd = n2 / (n1 + n2)
            b = lbd * b1
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = D(subsample1, grid, s) - D(subsample2, grid, s) # ngrid x N-b+1
            
            selected_D_sub = D_sub + recentering_function
            
            # take supremum over support
            return (b) ** (0.5) * (selected_D_sub).max(0)[:, None]
                
        elif resampling == 'bootstrap' :
            bootstrap = self.bootstrap

            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nsamp)
            btspsample2 = bootstrap(sample2, b2, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            selected_D_btsp = D_btsp + recentering_function
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (selected_D_btsp).max(0)[:, None]
        
        elif resampling == "paired_bootstrap" :
            
            paired_bootstrap = self.paired_bootstrap
            
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            selected_D_btsp = D_btsp + recentering_function
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (selected_D_btsp).max(0)[:, None]
    
        else :
            print("check resampling parameter")
      
        
    def plot_D(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        Parameters
        ==========
        save  : boolean
            if True, save figure
        titel : str
            name of the figure to save
        
        """
        
        
        sample1 = self.sample1
        sample2 = self.sample2
        
        grid    = self.grid
        s       = self.s
        D       = self.D

        D1 = D(sample1,grid,s)
        D2 = D(sample2,grid,s)
        
        # if s == 1 :
        #     ylabel = "Cumulative distribution function"
        # elif s > 1 :
        #     ylabel = str(s) + "-time cumulated distribution function"
            
        plt.figure()
        plt.plot(grid, D1, label = label1)
        plt.plot(grid, D2, label = label2, ls = 'dashed')
        plt.xlabel(xlabel)
        # plt.ylabel(ylabel)
        plt.legend()
        
        if save == True:
            plt.savefig(title)

        plt.show()
        
        
#########################################################
################# Hong and Li ###########################
#########################################################

        
class test_sd_NDM :
  
    """

    Stochastic Dominance Testing via Numerical Delta Method (Hong and Li, 2019)
   
    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0
    
    <=> 
    
    H0 : d_s <= 0 
    H1 : d_s > 0
    
    Parameters
    ----------
    sample1 : np.array (1-d)
    sample2 : np.array (1-d)
                input samples should be 1-dim np.array
    
    grid    : int
                # of grid points
    s       : int
                Order of Stochastic Dominance
    b1      : int
                resampling size of the sample1.
    b2      : int
                resampling size of the sample2.
    resampling : str
                resampling should be one of 'subsampling' and 'bootstrap'  
                
                'subsampling' -> use sumbsampling method suggested by LMW
                'bootstrap'   -> use recentered bootstrap method as in BD and LMW
                
    nsamp   : int
                # of bootstrap statistics for bootstrap distribution
    
    epsilon : float
                Tuning parameter for NDM (default value: r_N^(-1/16))
                
    form    : str
                Specifying functional forms as KS, L1, or L2
                
                'KS' : Kolmogorov-Smirnov type test statistics (Supremum type)
                'L1' : L1 type test statistics (Integral type)
                'L2' : L2 type test statistics (Integral type)

    """
    
    def __init__(self, sample1, sample2, ngrid, s, b1, b2, resampling, nsamp = 200, epsilon = None, form = "KS"):
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.b1          = b1
        self.b2          = b2
        self.nsamp       = nsamp
        
        n1 = sample1.shape[0] 
        n2 = sample2.shape[0] 
        r_N =  (n1 * n2 / (n1 + n2)) ** (0.5)
        
        self.n1          = sample1.shape[0]   
        self.n2          = sample2.shape[0]   
        self.resampling  = resampling

        # Tuning Parameters

        if epsilon == None:
            epsilon = r_N ** (-1/16);
        
        self.epsilon     = epsilon
        
        # KS-type statistics is default
        self.form        = form

        # setting grid
        minx = min(sample1.min(0), sample2.min(0))
        maxx = max(sample1.max(0), sample2.max(0))     
        grid = np.linspace(minx,maxx, num=ngrid)
        
        self.grid = grid
        
        
    def testing(self) :
        
        """
        
        Returns
        -------
            None
            In here, we print the result
            
        test_stat   : float
                the value of the test statistic
        test_stat_b : numpy array
                values of the resampled test statistics
        pval        : float
                p-value of the test
        
        
        
        """
        
        sample1     = self.sample1
        sample2     = self.sample2
        ngrid       = self.ngrid
        s           = self.s
        b1          = self.b1
        b2          = self.b2
        resampling  = self.resampling
        nsamp       = self.nsamp
        epsilon     = self.epsilon
        
        start_time = time.time()
                
        # grid
        grid  = self.grid
        
        # Estimation
        test_stat   = self.T_N(sample1, sample2, grid, s)
        self.test_stat = test_stat
        
        test_stat_b = self.resampled_stat(sample1, sample2, grid, s, b1, b2, resampling, nsamp)
        pval = (test_stat_b >= test_stat).mean(0) 
        
        if s == 1:
            Torder = 'First Order SD'
        elif s == 2:
            Torder = 'Second Order SD'
        elif s == 3:
            Torder = 'Third Order SD'
        else:
            Torder = str(s) + 'th order SD'
    
    

        print('\n#-------------------------------------------#')    
        print('Testing Stochastic Dominance by Numerical Delta Method (Hong and Li 2018)',
          '\n\n* H0 : sample1 ', Torder, 'sample2')
        
        print('* Resampling method \t:', resampling)
        print('\n')
        
        print('* # (sample1) \t\t = %6d' % sample1.shape[0],
          '\n* # (sample2) \t\t = %6d\n' % sample2.shape[0])
        print('* # ('+ resampling + '1) \t = %6d' % b1,
          '\n* # ('+ resampling + '2) \t = %6d\n' % b2)
        print('#-------------------------------------------#\n')    
        print('* SD order \t\t = %6d' % s,
          '\n* # of grid points \t = %6d\n' % ngrid)
        print('# Tuning parameter ---------')
        print('# epsilon = %6d' % epsilon)
        print('#-------------------------------------------#\n')    
        print('* Test Result *\n')    

        print('* Test statistic \t = %5.4f' % test_stat)
        print('* p-value \t\t = %5.4f\n' % pval)
        print('#-------------------------------------------#')    
        et = time.time() - start_time
        print('\n* Time elapsed : %5.2f Sec' % et)            
        
        self.result = {'test_stat'   : test_stat[0] ,
                   'test_stat_b' : np.squeeze(test_stat_b),
                   'pval'        : pval[0]
                   }
        
    def operator(self, sample, grid, s):
        
        # np.apply_along_axis(operator, axis, a, b)
        # apply operator (in here, op.lt) along axis
        # op.lt(a,b) <=> a < b
        # op.sub(a,b) <=> a - b
            
        return (np.apply_along_axis(op.lt, 1, sample, grid) * (- np.apply_along_axis(op.sub, 1, sample, grid)) ** (
            s - 1) / math.factorial(s - 1))
    
    def D(self, sample, grid, s):
        
        if sample.ndim != 3:       
            sample = sample[:,None,None]
        
        # We take average here.
        # if s = 1, D equals to ecdf
        
        operator = self.operator
        D = operator(sample, grid, s).mean(0)
        
        return D 
    
    def phi(self, theta):
    
        """ Set and calculate functional phi """
        
        form      = self.form
        
        # For original sample
        if theta.ndim == 1:
            if   form == 'KS':
                return (theta).max(0)
            elif form == 'L1':
                return np.trapz((theta>0) * theta, axis = 0)
            elif form == 'L2':
                return np.trapz(((theta>0) * theta) ** 2, axis = 0)
            else :
                print("Functional form should be KS or integral")
        
        # For resampled sample
        else:
            if   form == 'KS':
                return (theta).max(0)[:,None]
            elif form == 'L1':
                return np.trapz((theta>0) * theta, axis = 0)[:,None]
            elif form == 'L2':
                return np.trapz(((theta>0) * theta) ** 2, axis = 0)[:,None]
            else :
                print("Functional form should be 'KS', 'L1' or 'L2'")

    def T_N(self, sample1, sample2, grid, s):
            
        """
            
        test statsitics
        
        S_jf in BD
            
        """
        
        D    = self.D
        form = self.form
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = D(sample1, grid, s) - D(sample2, grid, s)
        self.D_s = D_s
        
        if form == "KS":
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)
        elif form == "L1":
            return (n1 * n2 / (n1 + n2)) ** (0.5) * np.trapz((D_s>0) * D_s, axis = 0)
        elif form == "L2":
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_s>0) * D_s)**(2), axis = 0)
                
        
    def subsampling(self, sample, subsize, nsub):
        
        """
        
        Subsampling
        
        Parameters
        ----------
        sample  : np.array (1-d)
        subsize : number of samples chosen by subsampling
        nsum    : number of subsample test statistics
        
        """
        
        
        # nsub : # of subsamples
        # subindex : subsample size x 1 x nsub
        subindex = (np.arange(0, subsize, 1)[:, None] + np.arange(0, nsub, 1)[:, None].T)[:, None]
    
        # subsample : subsample size x 1 x nsub x 1
        subsample = sample[subindex]
        return subsample
    
    def bootstrap(self, sample, btspsize, nbtsp):
        
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample = sample[btspindex]
        return btspsample
        
    
    def paired_bootstrap(self, sample1, sample2, btspsize, nbtsp):
        
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
        btspindex = np.array([np.random.randint(n, size=btspsize) for _ in np.arange(nbtsp)]).T[:,None,:]
        # btspsample : bootstrap size x 1 x nbtsp x 1 
        btspsample1 = sample1[btspindex]
        btspsample2 = sample2[btspindex]
        
        return btspsample1, btspsample2 
    
    def NDM(self, D_b_recentered):
    
        """
        
        Attaining aymptotic distributions by Numerical Delta Method (Hong and Li, 2018)
                
        Parameters
        ----------
            D_b_recentered : numpy array ()
            
        
        Results
        ----------
        the following variables are generated 
        
        phi_prime : numpy array (ngrid x 1)
        
        """
        
        D_s     = self.D_s
        n1      = self.n1
        n2      = self.n2
        epsilon = self.epsilon
        form    = self.form
        phi     = self.phi
        
        phi_D   = phi(D_s)
        
        # estimation
        r_N = ((n1*n2) / (n1 + n2))**(0.5)

        #  First order Approximation
        if   form == "KS" or form == "L1" :
            return (phi(D_s + epsilon*r_N*(D_b_recentered)) - phi_D) / epsilon
        
        #  Higher order Approximation
        elif form == "L2":
            return (phi(D_s + epsilon*r_N*(D_b_recentered)) - phi_D) / (epsilon**(2))
        
    
    
    def resampled_stat(self, sample1, sample2, grid, s, b1, b2, resampling, nsamp):
        
        D = self.D
        NDM = self.NDM
        
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
            
        if resampling == 'subsampling' :
            
            subsampling = self.subsampling

            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            lbd = n2 / (n1 + n2)
            b = lbd * b1
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = D(subsample1, grid, s) - D(subsample2, grid, s) # ngrid x N-b+1

            # Approximating asymptotic distribution via NDM                        
            return NDM(D_sub)
                
        elif resampling == 'bootstrap' :
            bootstrap = self.bootstrap

            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nsamp)
            btspsample2 = bootstrap(sample2, b2, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # Approximating asymptotic distribution via NDM
            return NDM(D_btsp)
        
        elif resampling == "paired_bootstrap" :
            
            paired_bootstrap = self.paired_bootstrap
            
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nsamp)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = D(btspsample1, grid, s) - D(sample1, grid, s) 
            D2_recentered = D(btspsample2, grid, s) - D(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered

            # Approximating asymptotic distribution via NDM            
            return NDM(D_btsp)
    
        else :
            print("check resampling parameter")
      
        
    def plot_D(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        Parameters
        ==========
        save  : boolean
            if True, save figure
        titel : str
            name of the figure to save
        label 1,2 : str
            the name of legend
        xlabel : str
            the label of x-axis.
        
        """
        
        
        sample1 = self.sample1
        sample2 = self.sample2
        
        grid    = self.grid
        s       = self.s
        D       = self.D

        D1 = D(sample1,grid,s)
        D2 = D(sample2,grid,s)
        
        # if s == 1 :
        #     ylabel = "Cumulative distribution function"
        # elif s > 1 :
        #     ylabel = str(s) + "-time cumulated distribution function"
            
        plt.figure()
        plt.plot(grid, D1, label = label1)
        plt.plot(grid, D2, label = label2, ls = 'dashed')
        plt.xlabel(xlabel)
        # plt.ylabel(ylabel)
        plt.legend()
        
        if save == True:
            plt.savefig(title)

        plt.show()           