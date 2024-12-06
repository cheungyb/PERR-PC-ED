** PERR with recurrent events
** Simulation Settings for scenario (A)
** Equivalence between PERR_Cox and PERR*_Cox
** 2 levels of baseline incidence; 3 levels of heterogeneity

clear
set more off

capture program drop PERR_AG_Program
program define PERR_AG_Program, rclass

version 15.1
 syntax [, n(integer 3000) beta0(real -7.0) beta(real 0.5) var(real 0.5) ka(real 1.25) kaC(real 1.25)   ///
			c0(real -8.0) cx(real 0.5) maxE(integer 100) bn(integer 2) ]
			
drop _all
local obs=`n'
set obs `obs'

local sdE=sqrt(`var')

gen id=_n
*time in days
gen tau=runiform(200,300)

*time constant covariates
gen z=rbinomial(1,0.5)
gen x=rbinomial(1,0.5)
gen exit=tau

gen v=rnormal(0,`sdE')  /* sd^2 = 0, 0.5, 1   */

*time to receive treatment
gen lamC = 1.0*exp(`c0'+`cx'*x+0.5*z+rnormal(0,0.1))
gen ct=(-ln(runiform())/lamC)^(1/`kaC')
replace ct=exit if ct>exit

*time to recurrent events
*true log HR
local beta1=log(`beta')
*hazard of events: Weibull distribution with ka>1
gen lam = 0.5*exp(`beta0'+`cx'*x+0.5*z+v)
gen t1=(-ln(runiform())/lam)^(1/`ka')
replace t1= (ct^`ka'-ln(runiform())/(exp(`beta1'*(ct<exit))*lam))^(1/`ka') if t1>ct 
replace t1=exit if t1>exit

qui{
forvalues i=2(1)`maxE'{
local j=`i'-1
gen t`i' = (t`j'^`ka'-ln(runiform())/(exp(`beta1'*(ct<exit)*(t`j'>ct))*lam))^(1/`ka')
replace t`i' = (ct^`ka'-ln(runiform())/(exp(`beta1'*(ct<exit))*lam))^(1/`ka')  if t`i'>ct & t`j'<ct

replace t`i'=exit if t`i'>exit

}
}
cap drop trt
gen trt=0
replace trt=1 if ct<exit
drop lam* v tau

reshape long t ,i(id) j(temp)
drop temp 

duplicates drop
drop if t==0

gen event=1
replace event=0 if  t==exit 

by id (t),sort: gen n=_n
by id (t),sort: gen N=_N

gen ep=1
replace ep=2 if n==N
expand ep
drop n N

sort id t ep
by id t ep,sort: gen n=_n
by id t ep,sort: replace t=ct if n==2

drop n ep
drop if t>exit
duplicates drop

gen start=0
gen end=t
by id (t), sort: replace start=end[_n-1] if _n>1
drop t

*time-varying variable 'treat'
gen treat=0
replace treat=1 if start>=ct

*observed time-constant covariate: 'x' 'z'

********************************************
** AG model 1, Time-on-study as timescale
tab event
stset end, failure(event) exit(t .) id(id) enter(start) 
stcox treat x z, hr cluster(id)

matrix table_Treat = r(table)
return scalar coef_Treat = table_Treat[1,1]
return scalar se_Treat = table_Treat[2,1]
return scalar pvalue_Treat = table_Treat[4,1]
return scalar cp_Treat = cond(`beta'>=table_Treat[5,1]&`beta'<=table_Treat[6,1],1,0)

stset, clear
drop treat 

********************************************

sort id start end
by id: gen n=_n

reshape wide start end event, i(id) j(n)

gen cut=ct
replace cut=exit if cut>=exit

gen entry0=0
stset cut, failure(trt) enter(entry0) 

gen inclusive= (trt==0)
qui sttocc_rev, match(z) n(1) nodots

gen index = ct
by _set (_case), sort: replace index=index[_N]

by _set, sort: gen Nset=_N
drop if Nset==1

drop Nset _set _time _case  cut
drop inclusive entry0  

reshape long start end event, i(id) j(temp)
drop temp
drop if event==.

expand 2 if trt==0

sort id start end
by id start: gen n=_n
replace end=index if n==1 & end>index & start<index & trt==0
replace start=index if n==2 & end>index & start<index & trt==0
drop n
duplicates drop

replace event=0 if end==index & trt==0

gen period=cond(end<=index,0,1)

sort id end
by id: gen n=_n

qui su ct if trt==1 & n==1,det
return scalar ttotrt_median=r(p50)

count if trt==1 & n==1
return scalar pop1=r(N)

qui su event if trt==0 & period==0
return scalar  neverE0 = r(sum)
qui su event if trt==1 & period==0
return scalar  everE0 = r(sum)

qui su event if trt==0 & period==1
return scalar  neverE1 = r(sum)
qui su event if trt==1 & period==1
return scalar everE1 = r(sum)
drop n

*PERR analysis from here:	
gen pt=trt*period

*PERR from AG model
stset end , failure(event) exit(t .) id(id) enter(start) 
stcox  pt trt period, hr cluster(id)

matrix table_AG = r(table)
return scalar coef_AG = table_AG[1,1]
return scalar se_AG = table_AG[2,1]
return scalar pvalue_AG = table_AG[4,1]
return scalar cp_AG = cond(`beta'>=table_AG[5,1]&`beta'<=table_AG[6,1],1,0)

				 
*Cox model to first event:
gen post_event= (period==1 & event==1)
sort id end
by id: gen post=sum(post_event)
gen first1= (post==1) & (post_event==1)
tab first1 trt if period==1

drop post_event 

gen prior_event= (period==0 & event==1)
sort id end
by id: gen prior=sum(prior_event)
gen first0= (prior==1) & (prior_event==1)
tab first0 trt if period==0

drop prior_event 

stset end if period==0, failure(first0)  id(id) enter(start) 
stcox trt , hr
matrix table_Cox0 = r(table)

gen st=_st
gen d=_d
gen t=_t
gen t0=_t0

stset end if period==1, failure(first1)  id(id) enter(start) 
stcox trt , hr
matrix table_Cox1 = r(table)
			 
replace _st=st if _st==0
replace _d=d if _d==.
replace _t=t if _t==.
replace _t0=t0 if _t0==.
drop st-t0

* PERR from Cox model
stcox  pt trt period, hr cluster(id)

matrix table_Cox = r(table)
return scalar coef_Cox = table_Cox[1,1]
return scalar se_Cox = table_Cox[2,1]
return scalar pvalue_Cox = table_Cox[4,1]
return scalar cp_Cox = cond(`beta'>=table_Cox[5,1]&`beta'<=table_Cox[6,1],1,0)
			 
* bootstrapping to get SE of ratio of log(HR) from Cox1 and Cox0
return scalar coef_CoxRR = table_Cox1[1,1]/table_Cox0[1,1]

bootstrap CoxRR=r(ratio), reps(`bn')  cluster(id) idcluster(newpairid) nowarn nodots:  ///
		  Cox_Ratio newpairid trt  start end period  post first1  prior first0
  
matrix table_cRR = r(table)
return scalar se_CoxRR = exp(table_cRR[2,1])
return scalar pvalue_CoxRR = table_cRR[4,1]

matrix table_PI=e(ci_percentile)
return scalar cp_CoxRR = cond(`beta'>=exp(table_PI[1,1])&`beta'<=exp(table_PI[2,1]),1,0)

eret clear

end

