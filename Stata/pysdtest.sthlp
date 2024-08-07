{smcl}
{* *! version 1.0  31 Jul 2024}{...}
{viewerjumpto "Syntax" "pysdtest##syntax"}{...}
{viewerjumpto "Description" "pysdtest##description"}{...}
{viewerjumpto "Options" "pysdtest##options"}{...}
{viewerjumpto "Results" "pysdtest##results"}{...}
{viewerjumpto "Examples" "pysdtest##examples"}{...}
{viewerjumpto "References" "pysdtest##references"}{...}
{title:Title}

{phang}
{bf:pysdtest} {hline 2} testing for stochastic dominance


{marker syntax}{...}
{title:Syntax}

{p 8 17 3}
{cmd:pysdtest} {varlist}(max=2) [if] [, by(varname) SWitch resampling(string) 
approach(string) Functional(string) ngrid(integer 100) s(integer 1) b1(integer 0) b2(integer 0) 
nboot(integer 200) c(real 0.75) a(real 0.1) eta(real 0.000001) epsilon(real 0) alpha(real 0.05)]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}

{synopt:{opt by}}Specify a binary variable for dividing the sample. If the option by( ) is not given, there must be two (continuous) variables for performing testing.  {p_end}

{synopt:{opt sw:itch}}Whether to switch an order of the sample when the option by( ) is specified.{p_end}

{synopt:{opt resampling}}Specify a resampling method. The method must be one of "bootstrap", "subsampling", or "paired_bootstrap". If the option is not given, "bootstrap" is used as a default option.  To use "subsampling", the number of subsample of first and second variables must be given by options b1( ) and b2( ), respectively. To use "paired_bootstrap", variables need to be arranged in a correct pair and have the same number of observations. {p_end}

{synopt:{opt approach}}Specify a testing method. If the option is not given, Kolmogorov-Smirnov type test statistic is used as in Barrett and Donald (2003) and Linton et al. (2005). If "contact" is given, the contact set approach by Linton et al. (2010) is performed. If "SR" is given, the selective recentering approach by Donald and Hsu (2016) is performed. If "NDM" is given, the numerical delta method is performed. {p_end}

{synopt:{opt f:unctional}}Specify a type of the functional for using the numerical delta method when "NDM" is given for the option approach( ). It must be one of "KS", "L1", and "L2". If the option is not given, "L1" is used as a default option. {p_end}

{synopt:{opt ngrid}}Set the number of grid points for calculating test statistics. 100 is set as a default. {p_end}

{synopt:{opt s}}Set the stochastic order of the null hypothesis. 1 is set as a default. It must be a positive integer. {p_end}

{synopt:{opt b1}}Set the number of subsample for the first variable. {p_end}

{synopt:{opt b2}}Set the number of subsample for the second variable. {p_end}

{synopt:{opt nboot}}Set the number of bootstrapping. {p_end}

{synopt:{opt c}}Set a value of the tuning parameter for the contact set approach. 0.75 is set as a default value. {p_end}

{synopt:{opt a}}Set a value of the tuning parameter for the selective recentering approach. 0.1 is set as a default value. {p_end}

{synopt:{opt eta}}Set a value of the tuning parameter for the selective recentering approach. 10e^{-6} is set as a default value. {p_end}

{synopt:{opt epsilon}}Set a value of the tuning parameter for the numerical delta method. 0.75 is set as a default value. {p_end}

{synopt:{opt alpha}}Set the significance level for the statistical test. 0.05 is set as a default value. {p_end}

{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd} {cmd: pysdtest} performs statistical tests for the null hypothesis of stochastic dominance (SD) that is the first sample (variable) s-th order stochastically dominates the second one. s can be set as an arbitrary stochastic order.

{pstd} Supported testing procedures include, but are not limited to, tests through the least favorable case approach of Barrett and Donald (2003), subsampling approach of Linton et al. (2005), contact set approach of Linton et al. (2010) and selective
recentering method of Donald and Hsu (2016). In addition, the numerical delta method (NDM; Dümbgen 1993; Fang and Santos 2019; Hong and Li 2018) can be used to compute the critical value. 

{pstd} {cmd:pysdtest} also supports the combination of various resampling methods, test statistics and procedures for approximating the limiting distribution.

{pstd} The command is based on the python package PySDTest. Therefore Stata with version higher than 16.0, Python 3 software and installation of the Python package are necessary for using the command.

{pstd} See {browse "https://arxiv.org/pdf/2307.10694":Lee and Whang (2024)} on arXiv for further reference on concepts/tests of SD and details/guidance about using the Stata/Python package.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:pysdtest} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N1)}}number of observations of the first variable {p_end}
{synopt:{cmd:r(N1)}}number of observations of the second variable {p_end}
{synopt:{cmd:r(b1)}}number of subsamples of the first variable {p_end}
{synopt:{cmd:r(b2)}}number of subsamples of the second variable {p_end}
{synopt:{cmd:r(s)}}stochastic order{p_end}
{synopt:{cmd:r(alpha)}}significance level{p_end}
{synopt:{cmd:r(ngrid)}}number of grid points{p_end}
{synopt:{cmd:r(test_stat)}}value of the test statistic{p_end}
{synopt:{cmd:r(p_val)}}p-value{p_end}
{synopt:{cmd:r(critical_val)}}critical-value{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(resampling)}}resampling method{p_end}
{synopt:{cmd:r(approach)}}testing approach{p_end}
{synopt:{cmd:r(functional)}}type of functional for the numerical delta method{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(grid)}}matrix of grid values{p_end}
{synopt:{cmd:r(limit_dist)}}matrix of the approximated limit distribution by resampling{p_end}


{marker examples}{...}
{title:Examples}

{phang}{cmd:. pysdtest var1 var2}{p_end}

{phang}{cmd:. pysdtest var1 var2, resampling("subsampling") b1(100) b2(200)}{p_end}

{phang}{cmd:. pysdtest var1 var2, resampling("paired_bootstrap") approach("contact") c(0.5)}{p_end}

{phang}{cmd:. pysdtest var, by(binary_var)}{p_end}

{phang}{cmd:. pysdtest var, by(binary_var) switch}{p_end}

{phang}{cmd:. pysdtest var, by(binary_var) approach("NDM") functional("KS") ngrid(100) nboot(500)}{p_end}

{marker references}{...}
{title: References}

{pstd} Barrett, Garry F., and Stephen G. Donald. "Consistent tests for stochastic dominance." Econometrica 71, no. 1 (2003): 71-104.

{pstd} Donald, Stephen G., and Yu-Chin Hsu. "Improving the power of tests of stochastic dominance." Econometric Reviews 35, no. 4 (2016): 553-585.

{pstd} Dümbgen, Lutz. "On nondifferentiable functions and the bootstrap." Probability Theory and Related Fields 95 (1993): 125-140.

{pstd} Fang, Zheng, and Andres Santos. "Inference on directionally differentiable functions." The Review of Economic Studies 86, no. 1 (2019): 377-412.

{pstd} Hong, Han, and Jessie Li. "The numerical delta method." Journal of Econometrics 206, no. 2 (2018): 379-394.

{pstd} Lee, Kyungho, and Yoon-Jae Whang. "PySDTest: a Python/Stata Package for Stochastic Dominance Tests." arXiv preprint arXiv:2307.10694 (2024).

{pstd}
Linton, Oliver, Esfandiar Maasoumi, and Yoon-Jae Whang. "Consistent testing for stochastic dominance under general sampling schemes." The Review of Economic Studies 72, no. 3 (2005): 735-765.

{pstd}
Linton, Oliver, Kyungchul Song, and Yoon-Jae Whang. "An improved bootstrap test of stochastic dominance." Journal of Econometrics 154, no. 2 (2010): 186-202.

