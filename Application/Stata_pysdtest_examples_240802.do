/*

Examples of using pysdtest

Authors: Kyungho Lee and Yoon-Jae Whang
Date: Aug 02. 2024.

*/

log using pysdtest_examples.smcl


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

* Run pysdtest 
pysdtest X1 X2

* Run pysdtest with subsampling 
pysdtest X1 X2, resampling("subsampling") b1(40) b2(40)

* Run pysdtest with paired bootstrapping and contact set approach
pysdtest X1 X2, resampling("paired_bootstrap") approach("contact") c(0.5)


******************
* Using auto data
******************
sysuse auto, clear

* Run pysdtest with by( ) option
pysdtest price, by(foreign)

* Run pysdtest with by( ) option and switch
pysdtest price, by(foreign) switch

* Run pysdtest with by( ) option, NDM
pysdtest price, by(foreign) approach("NDM") functional("KS") ngrid(200)
