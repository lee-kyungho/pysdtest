"""

@author: Kyungho Lee
Latest Update at July 30, 2024

Python Package for Stochastic Dominance Tests

"""

# Import modules
from pysdtest.Resampling_Functions import bootstrap, subsampling, paired_bootstrap
from pysdtest.Test_Stat_Functions import CDF, set_grid
import numpy as np
import matplotlib.pyplot as plt
import time

#%%
class test_sd :
    
    """

    Stochastic Dominance Test by Barett and Donald (2003 ECMA, hereafter BD), Linton, Maasoumi and Whang (2005 RESstud, hereafter LMW)

    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0
    
    Parameters
    ==============================
    sample1 : np.array (1-d)
    sample2 : np.array (1-d)
                input samples should be 1-dim np.array
    
    ngrid   : int
                # of grid points
    s       : int
                Order of Stochastic Dominance
    b1      : int
                resampling size of the sample1.
    b2      : int
                resampling size of the sample2.
    resampling : str
               
                resampling should be one of 'subsampling', 'bootstrap', and 'paired_bootstrap'   
                'subsampling'       -> use sumbsampling method as in LMW
                'bootstrap'         -> use recentered bootstrap method
                'paired_bootstrap'  -> use recentered bootstrap method by resampling a pair (X,Y) to allow dependency of paired observations.
                
    nboot   : int
                # of bootstrap statistics for the bootstrap distribution
    
    """
    
    def __init__(self, sample1, sample2, ngrid, s, resampling, b1 = None, b2 = None, nboot = 200, alpha = 0.05, quiet = False):
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.resampling  = resampling
        self.nboot       = nboot
        self.alpha       = alpha
        self.quiet       = quiet
        
        # set grid
        samples = [sample1, sample2]
        grid = set_grid(samples, ngrid)
        self.grid = grid
        
        if resampling == "bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
        
        elif resampling == "paired_bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
            
        self.b1          = b1
        self.b2          = b2

        
    def testing(self) :
        
        """
        
        Returns
        ==============================            
        test_stat   : float
                the value of the test statistic
        test_stat_b : numpy array
                values of the resampled test statistics
        critival_val: float
                critical value of the test
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
        nboot       = self.nboot
        quiet       = self.quiet
        
        start_time = time.time()
        
        # Estimation
        test_stat   = self.T_N()
        test_stat_b = self.resampled_stat()
        
        # calculate critical value and p-value
        alpha = self.alpha
        critival_val = np.quantile(test_stat_b, 1 - alpha)
        p_val = (test_stat_b >= test_stat).mean(0) 
        
        if s == 1:
            Torder = 'first order SD'
        elif s == 2:
            Torder = 'second order SD'
        elif s == 3:
            Torder = 'third order SD'
        else:
            Torder = str(s) + 'th order SD'

        et = time.time() - start_time
   
        if quiet == False:
            print('\n#--- Testing for Stochastic Dominance  -----#\n')    
            print('* H0 : sample1', Torder, 'sample2\n')
            print('#-------------------------------------------#\n')    
            print('*** Test Setting ***')    
            print('* Resampling method \t =', resampling)
            print('* SD order       \t = %6d' % s)
            print('* # of (sample1) \t = %6d' % sample1.shape[0],
            '\n* # of (sample2)   \t = %6d' % sample2.shape[0])
            if self.resampling == 'subsampling':
                print('* # of (subsample1) \t = %6d' % b1)
                print('* # of (subsample2) \t = %6d\n' % b2)
            else:
                print('* # of bootstrapping \t = %6d' % nboot)
                print('* # of grid points \t = %6d\n' % ngrid)
            print('#-------------------------------------------#\n')    
            print('*** Test Result ***')    
            print('* Test statistic \t = %5.4f' % test_stat)
            print('* Significance level \t = %5.2f' % alpha)
            print('* Critical-value \t = %5.4f' % critival_val)
            print('* P-value        \t = %5.4f' % p_val)
            print('* Time elapsed   \t = %5.2f Sec' % et)            
        
        # save the results
        self.result = {'test_stat'   : test_stat.reshape(-1)[0],
                   'test_stat_b' : np.squeeze(test_stat_b),
                   'critical_val' : critival_val,
                   'p_val'        : p_val[0]
                   }
        
    
    def T_N(self):
           
        sample1 = self.sample1
        sample2 = self.sample2
        grid    = self.grid
        s       = self.s
        
        """
            
        test statsitics
        
        d_s in LMW
        S_jf in BD
            
        """
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = CDF(sample1, grid, s) - CDF(sample2, grid, s)
        
        if D_s.ndim == 1:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)
        
        else:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)[:,None]
        
    def resampled_stat(self):
        
        sample1   = self.sample1
        sample2   = self.sample2
        grid      = self.grid
        s         = self.s
        b1        = self.b1
        b2        = self.b2
        resampling = self.resampling
        nboot     = self.nboot
        
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
            
        if resampling == 'subsampling' :
            
            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            lbd = n2 / (n1 + n2)
            b = lbd * b1
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = CDF(subsample1, grid, s) - CDF(subsample2, grid, s) # ngrid x N-b+1
            
            # take supremum over support
            return (b) ** (0.5) * (D_sub).max(0)[:, None]
                
        elif resampling == 'bootstrap' :
            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nboot)
            btspsample2 = bootstrap(sample2, b2, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (D_btsp).max(0)[:, None]
        
        elif resampling == "paired_bootstrap" :
                        
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (D_btsp).max(0)[:, None]
    
        else :
            print("check resampling parameter")
      
        
    def plot_CDF(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        
        Parameters
        ==============================
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

        D1 = CDF(sample1,grid,s)
        D2 = CDF(sample2,grid,s)
        
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
        
#%%        
class test_sd_contact:
  
    """

    Stochastic Dominance Testing by Linton, Song and Whang (2010 J.Econometrics, LSW)

    L-2 type Test Statistic 
    Contact-set estimation is used for constructing resampled statistics
    
    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0

    Parameters
    ==============================
    sample1 : np.array (1-d)
    sample2 : np.array (1-d)
                input samples should be 1-dim np.array
    
    ngrid    : int
                # of grid points
    s       : int
                Order of Stochastic Dominance
    b1      : int
                resampling size of the sample1.
    b2      : int
                resampling size of the sample2.
    resampling : str

                resampling should be one of 'subsampling', 'bootstrap', and 'paired_bootstrap'   
                'subsampling'       -> use sumbsampling method as in LMW
                'bootstrap'         -> use recentered bootstrap method
                'paired_bootstrap'  -> use recentered bootstrap method by resampling a pair (X,Y) to allow dependency of paired observations.
                
    nboot   : int
                # of bootstraping (default value: 200)
    
    c       : float
                Tuning parameter for contact-set estimation (default value: 0.75)
    
    """
    
    def __init__(self, sample1, sample2, ngrid, s, resampling, b1 = None, b2 = None, nboot = 200, c = 0.75, alpha = 0.05, quiet = False):
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.nboot       = nboot
        self.n1          = sample1.shape[0]   
        self.n2          = sample2.shape[0]   
        self.resampling  = resampling
        self.alpha       = alpha
        self.quiet       = quiet

        # setting grid
        samples = [sample1, sample2]
        grid = set_grid(samples, ngrid)
        self.grid = grid        
        
        if resampling == "bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
        
        elif resampling == "paired_bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
            
        self.b1          = b1
        self.b2          = b2
        
        # Tuning parameter
        self.c           = c
        N = (self.n1 + self.n2) / 2
        self.N           = N
        # sequence c_n for the contact set estimation
        self.seq_c       = c*np.log(np.log(N))/np.sqrt(N)
    
    def testing(self) :
        
        """
        
        Returns
        ==============================
            
        test_stat    : float
                the value of the test statistic
        test_stat_b  : numpy array
                values of the resampled test statistics
        critival_val : float
                critical value of the test
        p_val         : float
                p-value of the test
        
        """
        
        sample1     = self.sample1
        sample2     = self.sample2
        ngrid       = self.ngrid
        s           = self.s
        b1          = self.b1
        b2          = self.b2
        resampling  = self.resampling
        nboot       = self.nboot
        c           = self.c
        quiet       = self.quiet

        start_time = time.time()
        
        # Estimation
        test_stat   = self.T_N()
        test_stat_b = self.resampled_stat()
        
        # calculate critical value and p-value
        alpha = self.alpha
        critival_val = np.quantile(test_stat_b, 1 - alpha)
        p_val = (test_stat_b >= test_stat).mean(0) 
        
        
        if s == 1:
            Torder = 'first order SD'
        elif s == 2:
            Torder = 'second order SD'
        elif s == 3:
            Torder = 'third order SD'
        else:
            Torder = str(s) + 'th order SD'
    
        et = time.time() - start_time
        
        if quiet == False:
            print('\n#--- Testing for Stochastic Dominance  -----#\n')    
            print('* H0 : sample1', Torder, 'sample2')
            print('* Contact Set Approach\n')
            print('#-------------------------------------------#\n')    
            print('*** Test Setting ***')    
            print('* Resampling method \t =', resampling)
            print('* SD order       \t = %6d' % s)
            print('* # of (sample1) \t = %6d' % sample1.shape[0],
            '\n* # of (sample2)   \t = %6d' % sample2.shape[0])
            if self.resampling == 'subsampling':
                print('* # of (subsample1) \t = %6d' % b1)
                print('* # of (subsample2) \t = %6d\n' % b2)
            else:
                print('* # of bootstrapping \t = %6d' % nboot)
                print('* # of grid points \t = %6d\n' % ngrid)
            print("# Tuning parameter -------")
            print("* c              \t = %5.4f\n" % c)
            print('#-------------------------------------------#\n')    
            print('*** Test Result ***')    
            print('* Test statistic \t = %5.4f' % test_stat)
            print('* Significance level \t = %5.2f' % alpha)
            print('* Critical-value \t = %5.4f' % critival_val)
            print('* P-value        \t = %5.4f' % p_val)
            print('* Time elapsed   \t = %5.2f Sec' % et)            
    
        self.result = {'test_stat': test_stat.reshape(-1)[0] ,
                   'test_stat_b'  : np.squeeze(test_stat_b),
                   'critical_val' : critival_val,
                   'p_val'        : p_val[0]
                   }
        
        
    def T_N(self):
            
        sample1 = self.sample1
        sample2 = self.sample2
        grid    = self.grid
        s       = self.s
        
        """
            
        test statsitics
        
        d_s in LMW
        S_jf in BD
            
        """
                
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = CDF(sample1, grid, s) - CDF(sample2, grid, s)
        
        self.D_s  = D_s
        
        
        if D_s.ndim == 1:
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_s>0) * D_s)**(2),axis=0)
        
        else:
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_s>0) * D_s)**(2),axis=0)[:,None]
        

    def contact_set_estimation(self):
    
        """
        
        Contact-set estimation by Linton, Song and Whang (2010)
        This estimation is condcuted by using original samples.
        
        Parameters
        ==============================
            None
            
        
        Results
        ==============================        
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
        
    
    def resampled_stat(self):
        
        sample1    = self.sample1
        sample2    = self.sample2
        grid       = self.grid
        s          = self.s
        b1         = self.b1
        b2         = self.b2
        resampling = self.resampling
        nboot      = self.nboot
        
        contact_set_estimation = self.contact_set_estimation
        n1                     = self.n1
        n2                     = self.n2
        ngrid                  = self.ngrid 
        
        
        # contact-set estimation using the whole sample
        contact_set_estimation()        
        contact_set = self.contact_set
        ct_num      = self.ct_num
        
        
        if resampling == 'subsampling' :
            
            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            lbd = n2 / (n1 + n2)
            b = lbd * b1
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = CDF(subsample1, grid, s) - CDF(subsample2, grid, s) # ngrid x N-b+1
            
            # generate contact_set_rep for bootstrap sample
            # contact_set_rep is ngrid x nboot dimension
            
            contact_set_rep = np.repeat(contact_set,nsub) # (ngrid * nsub) x 1
            contact_set_rep = np.reshape(contact_set_rep, [ngrid, nsub]) # ngrid x nsub

            D_sub_contact = D_sub[contact_set_rep] # (ct_num * nsub) x 1
            D_sub_contact = np.reshape(D_sub_contact,[ct_num,nsub]) # ct_num x nsub
            
            
            # take integral on contact-set
            return (b) * np.trapz(((D_sub_contact>0) * D_sub_contact)**(2),axis=0)[:, None]
                
        elif resampling == 'bootstrap' :
            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nboot)
            btspsample2 = bootstrap(sample2, b2, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
        
            # generate contact_set_rep for bootstrap sample
            # contact_set_rep is ngrid x nboot dimension
            
            contact_set_rep = np.repeat(contact_set,nboot) # (ngrid * nboot) x 1
            contact_set_rep = np.reshape(contact_set_rep, [ngrid, nboot]) # ngrid x nboot

            D_btsp_contact = D_btsp[contact_set_rep]
            D_btsp_contact = np.reshape(D_btsp_contact,[ct_num,nboot])
            
            # take supremum over support
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_btsp_contact>0) * D_btsp_contact)**(2),axis=0)[:,None]
        
        elif resampling == "paired_bootstrap" :
                        
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # generate contact_set_rep for bootstrap sample
            # contact_set_rep is ngrid x nboot dimension
            
            contact_set_rep = np.repeat(contact_set,nboot) # (ngrid * nsamp) x 1
            contact_set_rep = np.reshape(contact_set_rep, [ngrid, nboot]) # ngrid x nsamp

            D_btsp_contact = D_btsp[contact_set_rep]
            D_btsp_contact = np.reshape(D_btsp_contact,[ct_num,nboot])
                        
            # take supremum over support
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_btsp_contact>0) * D_btsp_contact)**(2),axis=0)[:,None]
    
        else :
            print("check resampling parameter")
      
        
    def plot_CDF(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        
        Parameters
        ==============================
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

        D1 = CDF(sample1,grid,s)
        D2 = CDF(sample2,grid,s)
        
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
        
#%%                
class test_sd_SR :
  
    """

    Stochastic Dominance Testing by Donald and Hsu (2016, Econemtrics Review)

    KS-type Test Statistic with the Selective recentering approach
    
    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0
    
    Parameters
    ==============================
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
                'subsampling'       -> use sumbsampling method as in LMW
                'bootstrap'         -> use recentered bootstrap method
                'paired_bootstrap'  -> use recentered bootstrap method by resampling a pair (X,Y) to allow dependency of paired observations.
                
    nboot   : int
                # of bootstrap statistics for bootstrap distribution   
    a       : float
                Tuning parameter for selective-recentering approach. (default value: 0.1)
    eta     : float
                Tuning parameter for selective-recentering approach. (default value: 10^(-6))

    """
    
    def __init__(self, sample1, sample2, ngrid, s, resampling, b1 = None, b2 = None, nboot = 200, a = 0.1, eta = 10**(-6), alpha = 0.05, quiet = False):
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.nboot       = nboot
        self.alpha       = alpha
        self.quiet       = quiet
        
        n1 = sample1.shape[0] 
        n2 = sample2.shape[0] 
        
        self.n1          = n1   
        self.n2          = n2   
        self.resampling  = resampling
        self.a           = a

        # set grid
        samples = [sample1, sample2]
        grid = set_grid(samples, ngrid)
        self.grid = grid
        
        seq_a = -a* np.sqrt(np.log (np.log(n1+n2)))
        # Tuning parameter
        self.seq_a           = seq_a
        self.eta             = eta
        
        if resampling == "bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
        
        elif resampling == "paired_bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
            
        self.b1          = b1
        self.b2          = b2
        
    def testing(self) :
        
        """
        
        Returns
        ==============================
            
        test_stat    : float
                the value of the test statistic
        test_stat_b  : numpy array
                values of the resampled test statistics
        critival_val : float
                critical value of the test
        p_val         : float
                p-value of the test
        
        """
        
        sample1     = self.sample1
        sample2     = self.sample2
        ngrid       = self.ngrid
        s           = self.s
        b1          = self.b1
        b2          = self.b2
        resampling  = self.resampling
        nboot       = self.nboot
        a           = self.a
        eta         = self.eta
        quiet       = self.quiet
        
        
        start_time = time.time()
                
        # Estimation
        test_stat   = self.T_N()
        self.test_stat = test_stat
        
        test_stat_b = self.resampled_stat()
        test_stat_b[test_stat_b <= 10**(-6)] = 10**(-6)
        
        # calculate critical value and p-value
        alpha = self.alpha
        critival_val = np.quantile(test_stat_b, 1 - alpha)
        p_val = (test_stat_b >= test_stat).mean(0) 
        
        
        if s == 1:
            Torder = 'first order SD'
        elif s == 2:
            Torder = 'second order SD'
        elif s == 3:
            Torder = 'third order SD'
        else:
            Torder = str(s) + 'th order SD'
    
        et = time.time() - start_time

        if quiet == False:
            print('\n#--- Testing for Stochastic Dominance  -----#\n')    
            print('* H0 : sample1', Torder, 'sample2')
            print('* Selective Recentering Approach\n')
            print('#-------------------------------------------#\n')    
            print('*** Test Setting ***')    
            print('* Resampling method \t = ', resampling)
            print('* SD order       \t = %6d' % s)
            print('* # of (sample1) \t = %6d' % sample1.shape[0],
                '\n* # of (sample2)   \t = %6d' % sample2.shape[0])
            if self.resampling == 'subsampling':
                print('* # of (subsample1) \t = %6d' % b1)
                print('* # of (subsample2) \t = %6d\n' % b2)
            else:
                print('* # of bootstrapping \t = %6d' % nboot)
                print('* # of grid points \t = %6d\n' % ngrid)
            print('# Tuning paremeters -------------')
            print('* a              \t = %5.4f' % a)
            print('* eta            \t = %5.6f' % eta)        
            print('#-------------------------------------------#\n')    
            print('*** Test Result ***')    
            print('* Test statistic \t = %5.4f' % test_stat)
            print('* Significance level \t = %5.2f' % alpha)
            print('* Critical-value \t = %5.4f' % critival_val)
            print('* P-value        \t = %5.4f' % p_val)
            print('* Time elapsed   \t = %5.2f Sec' % et)            
    
        self.result = {'test_stat': test_stat.reshape(-1)[0],
                   'test_stat_b'  : np.squeeze(test_stat_b),
                   'critical_val' : critival_val,
                   'p_val'        : p_val[0]
                   }
        
    def T_N(self):
        
        """
            
        test statsitics
        
        S_jf in BD
            
        """
        
        sample1 = self.sample1
        sample2 = self.sample2
        grid    = self.grid
        s       = self.s        

        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = CDF(sample1, grid, s) - CDF(sample2, grid, s)
        self.D_s = D_s
        
        
        if D_s.ndim == 1:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)
        
        else:
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)[:,None]

    def selective_recentering(self):
    
        """
        
        Selective Recentering by Donald and Hsu (2014)
        This estimation is condcuted by using original samples.
        
        Parameters
        ==============================
        None
            
        
        Results
        ==============================        
        selected_set : numpy array (ngrid x 1)
        sl_num      : int 
        
        """
        
        seq_a = self.seq_a
        D_s   = self.D_s
        n1    = self.n1
        n2    = self.n2
        
        # estimation
        # pick indices for the 'contact set'
        selected_set = (np.sqrt(n2) * D_s) < seq_a
        
        # to identify whether the contact set is empty
        sl_num      = np.sum(selected_set)
        
        self.sl_num = sl_num        
        self.selected_set = selected_set
        
        return (n1 * n2 /(n1 + n2)) ** (0.5) * (D_s) * selected_set
        
    
    
    def resampled_stat(self):

        sample1    = self.sample1
        sample2    = self.sample2
        grid       = self.grid
        s          = self.s
        b1         = self.b1
        b2         = self.b2
        resampling = self.resampling
        nboot      = self.nboot
                


        selective_recentering = self.selective_recentering
        
        recentering_function = selective_recentering()
        
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
            
        if resampling == 'subsampling' :
            
            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            lbd = n2 / (n1 + n2)
            b = lbd * b1
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = CDF(subsample1, grid, s) - CDF(subsample2, grid, s) # ngrid x N-b+1
            
            selected_D_sub = D_sub + recentering_function
            
            # take supremum over support
            return (b) ** (0.5) * (selected_D_sub).max(0)[:, None]
                
        elif resampling == 'bootstrap' :
            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nboot)
            btspsample2 = bootstrap(sample2, b2, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            selected_D_btsp = D_btsp + recentering_function
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (selected_D_btsp).max(0)[:, None]
        
        elif resampling == "paired_bootstrap" :
                        
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            selected_D_btsp = D_btsp + recentering_function
            
            # take supremum over support
            return (n1 * n2 /(n1 + n2)) ** (0.5) * (selected_D_btsp).max(0)[:, None]
    
        else :
            print("check resampling parameter")
      
        
    def plot_CDF(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        Parameters
        ==============================
        save  : boolean
            if True, save figure
        titel : str
            name of the figure to save
        
        """
        
        
        sample1 = self.sample1
        sample2 = self.sample2
        
        grid    = self.grid
        s       = self.s

        D1 = CDF(sample1,grid,s)
        D2 = CDF(sample2,grid,s)
        
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
        
#%%
        
class test_sd_NDM :
  
    """

    Stochastic Dominance Testing via Numerical Delta Method (Dumbgen 1993, Fang and Santos 2019, Hong and Li, 2018)
   
    H0 : the sample1 's'th-order stochastically dominates the sample2
    H1 : the negation of H0
    
    Parameters
    ==============================
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
                'subsampling'       -> use sumbsampling method as in LMW
                'bootstrap'         -> use recentered bootstrap method
                'paired_bootstrap'  -> use recentered bootstrap method by resampling a pair (X,Y) to allow dependency of paired observations.
                
    nboot   : int
                # of bootstrap statistics for bootstrap distribution
    
    epsilon : float
                Tuning parameter for NDM (default value: r_N^(-1/16))
                
    form    : str
                Specifying functional forms as KS, L1, or L2
                
                'KS' : Kolmogorov-Smirnov type test statistics (Supremum type)
                'L1' : L1 type test statistics (Integral type)
                'L2' : L2 type test statistics (Integral type)

    """
    
    def __init__(self, sample1, sample2, ngrid, s, resampling, b1 = None, b2 = None, nboot = 200, epsilon = None, form = "L1", alpha = 0.05, quiet = False):
       
        self.sample1     = sample1
        self.sample2     = sample2
        self.ngrid       = ngrid
        self.s           = s
        self.nboot       = nboot
        self.alpha       = alpha
        self.quiet       = quiet

        n1 = sample1.shape[0] 
        n2 = sample2.shape[0] 
        r_N =  (n1 * n2 / (n1 + n2)) ** (0.5)
        
        self.n1          = n1   
        self.n2          = n2   
        self.resampling  = resampling

        # Tuning Parameters

        if epsilon == None:
            epsilon = r_N ** (-1/16)
        
        self.epsilon     = epsilon
        
        # KS-type statistics is default
        self.form        = form

        # set grid
        samples = [sample1, sample2]
        grid = set_grid(samples, ngrid)
        self.grid = grid
        
        
        if resampling == "bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
        elif resampling == "paired_bootstrap":        
            b1          = sample1.shape[0]
            b2          = sample2.shape[0]
            
        self.b1          = b1
        self.b2          = b2
        
        
    def testing(self) :
        
        """
        
        Returns
        ==============================
            
        test_stat    : float
                the value of the test statistic
        test_stat_b  : numpy array
                values of the resampled test statistics
        critival_val : float
                critical value of the test
        p_val         : float
                p-value of the test
        
        """
        
        sample1     = self.sample1
        sample2     = self.sample2
        ngrid       = self.ngrid
        s           = self.s
        b1          = self.b1
        b2          = self.b2
        resampling  = self.resampling
        nboot       = self.nboot
        epsilon     = self.epsilon
        quiet       = self.quiet
        form        = self.form
        
        start_time = time.time()
                        
        # Estimation
        test_stat   = self.T_N()
        self.test_stat = test_stat
        
        test_stat_b = self.resampled_stat()

        # calculate critical value and p-value
        alpha = self.alpha
        critival_val = np.quantile(test_stat_b, 1 - alpha)
        p_val = (test_stat_b >= test_stat).mean(0) 
                
        if s == 1:
            Torder = 'first order SD'
        elif s == 2:
            Torder = 'second order SD'
        elif s == 3:
            Torder = 'third order SD'
        else:
            Torder = str(s) + 'th order SD'

        et = time.time() - start_time

        if quiet == False:
            print('\n#--- Testing for Stochastic Dominance  -----#\n')    
            print('* H0 : sample1', Torder, 'sample2')
            print('* Numerical Delta Method')
            print('* Type of the Test Statistic \t = ', form)
            print('\n#-------------------------------------------#\n')    
            print('*** Test Setting ***')    
            print('* Resampling method \t =', resampling)
            print('* SD order       \t = %6d' % s)
            print('* # of (sample1) \t = %6d' % sample1.shape[0],
            '\n* # of (sample2)   \t = %6d' % sample2.shape[0])
            if self.resampling == 'subsampling':
                print('* # of (subsample1) \t = %6d' % b1)
                print('* # of (subsample2) \t = %6d\n' % b2)
            else:
                print('* # of bootstrapping \t = %6d' % nboot)
                print('* # of grid points \t = %6d\n' % ngrid)
            print('# Tuning paremeter -------------')
            print('* epsilon        \t = %5.4f\n' % epsilon)
            print('#-------------------------------------------#\n')    
            print('*** Test Result ***')    
            print('* Test statistic \t = %5.4f' % test_stat)
            print('* Significance level \t = %5.2f' % alpha)
            print('* Critical-value \t = %5.4f' % critival_val)
            print('* P-value        \t = %5.4f' % p_val)
            print('* Time elapsed   \t = %5.2f Sec' % et)            
    
        self.result = {'test_stat'   : test_stat.reshape(-1)[0],
                   'test_stat_b' : np.squeeze(test_stat_b),
                   'critical_val' : critival_val,
                   'p_val'        : p_val[0]
                   }        

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
                print("Functional form should be 'KS', 'L1' or 'L2'")
        
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

    def T_N(self):
        
        
        sample1 = self.sample1
        sample2 = self.sample2
        grid    = self.grid
        s       = self.s
        
        """
            
        test statsitics
        
        S_jf in BD
            
        """
        
        form = self.form
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
        
        D_s = CDF(sample1, grid, s) - CDF(sample2, grid, s)
        self.D_s = D_s
        
        if form == "KS":
            return (n1 * n2 / (n1 + n2)) ** (0.5) * (D_s).max(0)
        elif form == "L1":
            return (n1 * n2 / (n1 + n2)) ** (0.5) * np.trapz((D_s>0) * D_s, axis = 0)
        elif form == "L2":
            return (n1 * n2 / (n1 + n2)) * np.trapz(((D_s>0) * D_s)**(2), axis = 0)
                
            
    def NDM(self, D_b_recentered):
    
        """
        
        Attaining aymptotic distributions by Numerical Delta Method (Hong and Li, 2018)
                
        Parameters
        ==============================
        D_b_recentered : numpy array (ngrid x 1 x nbtsp (or nsub))
            
        
        Results
        ==============================        
        phi_prime : numpy array (ngrid x 1)
        
        """
        
        D_s     = self.D_s
        n1      = self.n1
        n2      = self.n2
        b1      = self.b1
        b2      = self.b2
        epsilon = self.epsilon
        form    = self.form
        phi     = self.phi
        resampling = self.resampling
        
        phi_D   = phi(D_s)
        
        # estimation
        
        if resampling == "subsampling":
            r_N = ((b1*b2) / (b1 + b2))**(0.5)            
            
        else:
            r_N = ((n1*n2) / (n1 + n2))**(0.5)

        #  First order Approximation
        if   form == "KS" or form == "L1" :
            return (phi(D_s + epsilon*r_N*(D_b_recentered)) - phi_D) / epsilon
        
        #  Higher order Approximation
        elif form == "L2":
            return (phi(D_s + epsilon*r_N*(D_b_recentered)) - phi_D) / (epsilon**(2))
        
    
    
    def resampled_stat(self):
        
        sample1    = self.sample1
        sample2    = self.sample2
        grid       = self.grid
        s          = self.s
        b1         = self.b1
        b2         = self.b2
        resampling = self.resampling
        nboot      = self.nboot
        NDM = self.NDM
        
        
        n1 = sample1.shape[0]
        n2 = sample2.shape[0]
            
        if resampling == 'subsampling' :
            
            # resampling via 'subsampling'
            
            # This part is to take into account different sample size (Linton, 2005)
            nsub = np.min([n1 - b1 + 1, n2 - b2 + 1])            
        
            subsample1 = subsampling(sample1, b1, nsub) # n1 x 1 x N-b+1
            subsample2 = subsampling(sample2, b2, nsub) # n2 x 1 x N-b+1
            
            D_sub = CDF(subsample1, grid, s) - CDF(subsample2, grid, s) # ngrid x N-b+1

            # Approximating asymptotic distribution via NDM                        
            return NDM(D_sub)
                
        elif resampling == 'bootstrap' :
            # resampling via 'boostrap'
                
            btspsample1 = bootstrap(sample1, b1, nboot)
            btspsample2 = bootstrap(sample2, b2, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered
            
            # Approximating asymptotic distribution via NDM
            return NDM(D_btsp)
        
        elif resampling == "paired_bootstrap" :
            
            paired_bootstrap = self.paired_bootstrap
            
            btspsample1, btspsample2 = paired_bootstrap(sample1,sample2, b1, nboot)
            
            if sample1.ndim != 3:       
                sample1 = sample1[:,None,None]
            if sample2.ndim != 3:             
                sample2 = sample2[:,None,None]
                    
            # recentering
            D1_recentered = CDF(btspsample1, grid, s) - CDF(sample1, grid, s) 
            D2_recentered = CDF(btspsample2, grid, s) - CDF(sample2, grid, s) 
            D_btsp = D1_recentered - D2_recentered

            # Approximating asymptotic distribution via NDM            
            return NDM(D_btsp)
    
        else :
            print("check resampling parameter")
      
        
    def plot_CDF(self, save = False, title = None, label1 = "sample1", label2 = "sample2", xlabel = "x"):
        
        """
        Parameters
        ==============================
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

        D1 = CDF(sample1,grid,s)
        D2 = CDF(sample2,grid,s)
        
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
