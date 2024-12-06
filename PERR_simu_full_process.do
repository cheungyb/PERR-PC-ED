** Simulation for PERR study
** 20 Nov 2024

clear 
set more off
version 15.1

cd "C:\Users\gmsmxm\NUS Dropbox\Xiangmei Ma\PERR2024\Recurrent\Simulation2\Code For Submission"

global reps=1000
run sttocc_rev.ado

********************************************************************************
** Simulation Settings for Table 1: Esimates from PERR, PERR_Cox and PERR_AG
** 2 levels of baseline incidence; 3 levels of heterogeneity
 
* 6 scenarios:
* b0 = -6.5 (high IR), -7.5 (low IR)
* HR = beta = 0.5
* var = 0, 0.5, 1

run PERR_AG_Program.ado
run coxRatio.ado
local hr = 0.5

foreach b0 in -6.5 -7.5 {
	foreach d0 in 0 0.5 1 {

simulate,reps($reps) seed(20240930) dots(10): ///
PERR_AG_Program, n(3000) beta0(`b0') var(`d0') beta(`hr') ka(1.25) bn(199)

foreach m in CoxRR Cox AG {
	qui su coef_`m'
	di "est_`m' = " r(mean)
	qui su cp_`m'
	di "CP_`m' = " r(mean)*100
	gen MSE_`m'=(exp(coef_`m')-`hr')^2
	qui su MSE_`m'
	di "RMSE_`m' = " sqrt(r(mean))
	di " "
		}	
	}
}

********************************************************************************
run PERR_EDT_program.ado

** Simulation Settings for Table 2, 3 and S1:
local hr = 0.5

foreach delta in  20 30 {
	foreach theta in 0.25 0.5 2 4 1 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta') beta0(-7) c0(-8) cx2(0.5) ka(1.25) 

if (`theta'<=1) {	
	foreach m in AG  AG_adj_all1  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last1  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	}
	
tab M
su coef_b3* if M==`delta'/10
su coef_b3* if abs(M-`delta'/10)<=1
su coef_b3*

	}
}

*************************
local delta = 30 
*** theta = 1/4 for 15 days followed by theta = 1/2 for another 15 days
simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(0.25) a0(0.5) 

foreach m in AG  AG_adj_all1  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
		
tab M
su coef_b3* if M==`delta'/10
su coef_b3* if abs(M-`delta'/10)<=1
su coef_b3*

*** theta = 4 for 15 days followed by theta = 2 for another 15 days
simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(4) a0(2) 

foreach m in AG  AG_adj_last1  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
		
tab M
su coef_b3* if M==`delta'/10
su coef_b3* if abs(M-`delta'/10)<=1
su coef_b3*


********************************************************************************
** Simulation Settings for Table S2 and S3:
local delta = 30

foreach hr in  1 2 {
	foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta')  

if (`theta'<=1) {	
	foreach m in AG  AG_adj_all1  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last1  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	}
	
    }
}

********************************************************************************
** For Table S4-S11:

local delta = 30
local hr = 0.5

** Simulation Settings for Table S4:

foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) sx(1) beta(`hr') gap(`delta') a1(`theta') a0(`theta')

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	}	
}

** Simulation Settings for Table S11:
foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) s(1) beta(`hr') gap(`delta') a1(`theta') a0(`theta')  

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	}
}

********************************************************************************
** Simulation Settings for Table S5:
foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta') beta0(-6.5) cx2(-0.5) 

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	}
}

** Simulation Settings for Table S6:
foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta')  beta0(-7.5) c0(-8.5) cx1(1) cx2(1)

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}	
	}	
}

********************************************************************************
** Simulation Settings for Table S7:
foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta') ka(1.0) beta0(-5.5)

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}	
	}	
}

** Simulation Settings for Table S8:
foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta') ka(1.5) beta0(-8.5)

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}	
	}	
}

** Simulation Settings for Table S9:
foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta') kaC(1.0) c0(-6.5)

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}	
	}
}

** Simulation Settings for Table S10:
foreach theta in 0.25 0.5 2 4 {

simulate,reps($reps) seed(20240930) dots(10): ///
	PERR_EDT_Program, n(3000) beta(`hr') gap(`delta') a1(`theta') a0(`theta') kaC(1.5) c0(-9.5)

if (`theta'<=1) {
	foreach m in AG  AG_adj_all2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}
	} 
else {
	foreach m in AG  AG_adj_last2 {
		qui su coef_`m'
		di "est_`m' = " r(mean)
		qui su cp_`m'
		di "CP_`m' = " r(mean)*100
		gen MSE_`m'=(exp(coef_`m')-`hr')^2
		qui su MSE_`m'
		di "RMSE_`m' = " sqrt(r(mean))
		di " "
		}	
	}
}
*End!
********************************************************************************

			
			
			
			
			
			
	