/*

Examples of using pysdtest

Authors: Kyungho Lee and Yoon-Jae Whang
Date: Aug 02. 2024.

*/

cd "../Application"

sjlog using pysdtest_examples_simul_data.smcl, replace

* Simulate the data
* Clear the existing dataset
clear

* Set the number of observations
set obs 100

* Set the seed
set seed 2024

* Generate a normally distributed variable
gen X1 = rnormal(0, 1)
gen X2 = rnormal(0, 1)

sjlog close, replace
sjlog using pysdtest_examples1.smcl, replace

* Run pysdtest 
pysdtest X1 X2

sjlog close, replace
sjlog using pysdtest_examples2.smcl, replace

* Run pysdtest with subsampling 
pysdtest X1 X2, resampling("subsampling") b1(40) b2(40)

sjlog close, replace
sjlog using pysdtest_examples3.smcl, replace

* Run pysdtest with paired bootstrapping and contact set approach
pysdtest X1 X2, resampling("paired_bootstrap") approach("contact") c(0.5)

sjlog close, replace
*************************************
* Using auto data with by( ) option *
*************************************
sjlog using pysdtest_examples4.smcl, replace

sysuse auto, clear

gen     foreign_str = "Domestic" if foreign == 0
replace foreign_str = "Foreign"  if foreign == 1

* Run pysdtest with by( ) option
pysdtest price, by(foreign_str)

sjlog close, replace
sjlog using pysdtest_examples5.smcl, replace

* Run pysdtest with by( ) option and switch
pysdtest price, by(foreign_str) switch

sjlog close, replace
sjlog using pysdtest_examples6.smcl, replace

* Run pysdtest with by( ) option, NDM
pysdtest price, by(foreign_str) switch ///
approach("NDM") functional("KS") ngrid(200)

sjlog close, replace
