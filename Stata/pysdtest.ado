program pysdtest, rclass
    version 16.0
    display "Running PySDTest"
    syntax varlist(max=2) [if] [, ///
	    by(varlist) ///
		SWitch ///
        resampling(string) ///
		approach(string) ///
		Functional(string) ///
		ngrid(integer 100) ///
        s(integer 1) ///
		b1(integer 0) ///
        b2(integer 0) ///
        nboot(integer 200) ///
		c(real 0.75) ///
		a(real 0.1) /// 
		eta(real 0.000001) ///
		epsilon(real 0) ///
		alpha(real 0.05)]
    
    *display "Syntax parsed successfully"
    	
    marksample touse
	
    // Validate resampling method
    if !inlist("`resampling'",  "", "bootstrap", "subsampling", "paired_bootstrap") {
        display as error "resampling must be one of: bootstrap, subsampling, or paired bootstrap"
         exit 198
     }

    // Validate subsampling sizes if method is subsampling
    if "`resampling'" == "subsampling" & (`b1' == 0 | `b2' == 0) {
        display as error "b1 and b2 must be specified for subsampling method"
        exit 198
    }

	local byvar `by'
    if "`byvar'" == "" {
		local var1: word 1 of `varlist'
		local var2: word 2 of `varlist'
		
		display "Sample1: `var1'"
		display "Sample2: `var2'"

		// Call Python function
		if ("`approach'" == "") | ("`approach'" == "contact") | ("`approach'" == "SR") | ("`approach'" == "NDM")  {
		python: run_test_sd("`var1'", "`var2'", "`touse'", "`approach'", `ngrid', `s', "`resampling'", `b1', `b2', `nboot', 								`alpha', `c', `a', `eta', `epsilon', "`functional'")
		}
		else {
			display as error "Specify a correct approach. It needs to be one of empty string '', 'contact', 'SR', or 'NDM'"
		}
		
    }

	* if by( ) is specified
	else {
	
	// Extract unique values of the by variable
	display "Groups:"
    levelsof `byvar' if `touse', local(levels)
    
    // Ensure there are exactly two levels
    local nlevels : word count `levels'
		if `nlevels' != 2 {
			display as error "by() variable must have exactly two levels"
			exit 198
		}
		
		local var: word 1 of `varlist'
		
		if "`switch'" != "" {
			local level1 : word 1 of `levels'
			local level2 : word 2 of `levels'
		}
		else {
			local level1 : word 2 of `levels'
			local level2 : word 1 of `levels'
		}

		display "Sample1: `var' of `level1'"
		display "Sample2: `var' of `level2'"

		// Call Python function
		if ("`approach'" == "") | ("`approach'" == "contact") | ("`approach'" == "SR") | ("`approach'" == "NDM")  {
		python: run_test_sd_by("`var'", "`byvar'", "`level1'", "`level2'", "`touse'", "`approach'", `ngrid', `s', 								"`resampling'", `b1', `b2', `nboot', `alpha', `c', `a', `eta', `epsilon', "`functional'")
		}
		else {
			display as error "Specify a correct approach. It needs to be one of empty string '', 'contact', 'SR', or 'NDM'"
		}	
	}
	
    // Return results
	// scalar
    return scalar N1  	       = r(N1)
    return scalar N2           = r(N2)
    return scalar b1  	       = r(b1)
    return scalar b2           = r(b2)
    return scalar s            = r(s)
    return scalar alpha        = r(alpha)
    return scalar ngrid        = r(ngrid)
    return scalar test_stat    = r(test_stat)
    return scalar p_val        = r(p_val)
    return scalar critic_val   = r(critic_val)

	// matrix
	return matrix grid = grid
	return matrix limit_dist = limit_dist
	
	// set global values
	global resampling        "`resampling'"
	global approach          "`approach'"
	global functional        "`functional'"

end

python:
import numpy as np
import pysdtest
from sfi import Data, Scalar, Matrix, Macro, SFIToolkit

def arrange_sample(var1, var2, touse):

	sample1 = np.array(Data.get(var1, missingval = np.nan))
	sample2 = np.array(Data.get(var2, missingval = np.nan))

	sample1 = sample1[~np.isnan(sample1)]
	sample2 = sample2[~np.isnan(sample2)]

	return sample1, sample2
	
def arrange_sample_by(var, byvar, level1, level2, touse):

	# Get all data
	all_data = np.array(Data.get([var, byvar, touse]))
	byvar_column = all_data[:,1]
	
	try:
		# try if the byvar is integer
		level1_filter = (byvar_column == int(level1))
		level2_filter = (byvar_column == int(level2))
	except:

		level1_filter = (byvar_column == level1)
		level2_filter = (byvar_column == level2)
		
	# Filter data for each level
	sample1 = all_data[level1_filter,0].astype(float)
	sample2 = all_data[level2_filter,0].astype(float)

	#print(f"Sample 1 size: {len(sample1)}")
	#print(f"Sample 2 size: {len(sample2)}")
	
	sample1 = sample1[~np.isnan(sample1)]
	sample2 = sample2[~np.isnan(sample2)]

	return sample1, sample2

def store_result(run_sd, result):

	N1 = run_sd.sample1.shape[0]
	N2 = run_sd.sample2.shape[0]
	b1 = run_sd.b1
	b2 = run_sd.b2
	s = run_sd.s
	alpha = run_sd.alpha
	grid  = run_sd.grid
	resampling = run_sd.resampling
	ngrid  = run_sd.ngrid
	limit_dist = result['test_stat_b']
	
	Scalar.setValue('r(N1)', N1)
	Scalar.setValue('r(N2)', N2)
	Scalar.setValue('r(b1)', b1)
	Scalar.setValue('r(b2)', b2)
	Scalar.setValue('r(alpha)', alpha)
	Scalar.setValue('r(ngrid)', ngrid)
	Scalar.setValue('r(s)', s)
	Scalar.setValue('r(test_stat)', result['test_stat'])
	Scalar.setValue('r(p_val)', result['p_val'])
	Scalar.setValue('r(critic_val)', result['critical_val'])
	
	# store matrices
	Matrix.store("grid", grid)
	Matrix.store("limit_dist", limit_dist)

	
def adjust_resampling(resampling, b1, b2):
	if resampling == "":
		resampling = "bootstrap"
	
	if resampling != "subsampling":
		b1 = None
		b2 = None
		
	return resampling, b1, b2

def run_test_sd(var1, var2, touse, approach, ngrid, s, resampling, b1, b2, nboot, alpha, c, a, eta, epsilon, functional):

	# arrange sample
	sample1, sample2 = arrange_sample(var1, var2, touse)

	# adjust arguments
	resampling, b1, b2 = adjust_resampling(resampling, b1, b2)

	# run sd test
	if approach == "":
		run_sd = pysdtest.test_sd(sample1, sample2, ngrid, s, resampling, b1, b2, nboot, alpha)
	
	elif approach == "contact":
		run_sd = pysdtest.test_sd_contact(sample1, sample2, ngrid, s, resampling,  b1, b2, nboot, c, alpha)
	
	elif approach == "SR":
		run_sd = pysdtest.test_sd_SR(sample1, sample2, ngrid, s, resampling, nboot, a, eta, alpha)
		
	elif approach == "NDM":
		if epsilon == 0:
			epsilon = None
			
		if functional == "":
			functional = "L1"
			
		run_sd = pysdtest.test_sd_NDM(sample1, sample2, ngrid, s, resampling, b1, b2, nboot, epsilon, functional, alpha)
	
	# perform testing
	run_sd.testing()
	result = run_sd.result

	# store the result
	store_result(run_sd, result)

# If byvar is specified

def run_test_sd_by(var, byvar, level1, level2, touse, approach, ngrid, s, resampling, b1, b2, nboot, alpha, c, a, eta, epsilon, functional):

	# arrange sample
	sample1, sample2 = arrange_sample_by(var, byvar, level1, level2, touse)

	# adjust arguments
	resampling, b1, b2 = adjust_resampling(resampling, b1, b2)

	# run sd test
	if approach == "":
		run_sd = pysdtest.test_sd(sample1, sample2, ngrid, s, resampling, b1, b2, nboot, alpha)
	elif approach == "contact":
		run_sd = pysdtest.test_sd_contact(sample1, sample2, ngrid, s, resampling,  b1, b2, nboot, c, alpha)
	elif approach == "SR":
		run_sd = pysdtest.test_sd_SR(sample1, sample2, ngrid, s, resampling, nboot, a, eta, alpha)
	elif approach == "NDM":
		if epsilon == 0:
			epsilon = None
		if functional == "":
			functional = "L1"
		run_sd = pysdtest.test_sd_NDM(sample1, sample2, ngrid, s, resampling, b1, b2, nboot, epsilon, functional, alpha)
	
	# perform testing
	run_sd.testing()
	result = run_sd.result

	# store the result
	store_result(run_sd, result)

end
