/*

Illustration of using pysdtest in Stata

Author: Kyungho Lee and Yoon-Jae Whang
Date: 240724

*/
cd "/Users/khlee/Dropbox/pysdtest/application"
import excel "BTC_daily_rr.xlsx", clear first
rename d_ln_Close BTC_daily_rr
save BTC_daily_rr, replace

import excel "sp_daily_rr.xlsx", clear first
rename d_ln_Close SP500_daily_rr
save sp_daily_rr, replace

merge using BTC_daily_rr
drop _merge

save bitcoin_sp500_daily_rr, replace

sjlog using bitcoin_sp500_daily_rr, replace
clear
use bitcoin_sp500_daily_rr, replace

* Specify the python environment
set python_exec /opt/anaconda3/bin/python
python query // Check the python version

* Run
pysdtest BTC_daily_rr SP500_daily_rr, ///
 resampling("subsampling") s(1) b1(1000) b2(900)
pysdtest BTC_daily_rr SP500_daily_rr, ///
 resampling("subsampling") s(2) b1(1000) b2(900)

pysdtest SP500_daily_rr BTC_daily_rr, ///
 resampling("subsampling") s(1) b1(900) b2(1000)
pysdtest SP500_daily_rr BTC_daily_rr, ///
 resampling("subsampling") s(2) b1(900) b2(1000)

sjlog close, replace

/*
* Run the stochastic dominance test on two variables
sysuse auto, clear

gen foreign_str = ""
replace foreign_str = "Domestic" if foreign == 0
replace foreign_str = "Foreign" if foreign == 1

discard
pysdtest mpg, by(foreign_str)
pysdtest price mpg, resampling("paired_bootstrap")  approach("contact") nboot(100)
help pysdtest


display r(test_stat)
display r(p_val)
display r(ngrid)
display r(b2)
display r(critic_val)
matrix list r(grid)
matrix list r(limit_dist)

display "$resampling"
display "$form"
display "$approach"

pysdtest price, by(foreign_str) resampling("subsampling") approach("NDM") s(1) nboot(200) form("KS") b1(20) b2(20)
*/
