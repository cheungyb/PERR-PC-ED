** PERR with recurrent events
** Simulation Settings for event dependent treatment time
** a0>1: events increase prob(trt==1); a0<1: events decrease prob(trt==1)
** two adjustments: 
** (1) shifting according to all prior events; (2) shifting according to last event

** change to be a1*h(trt) in first half gap, then a1*h(trt) in second half gap

clear
set more off

capture program drop PERR_EDT_Program
program define PERR_EDT_Program, rclass

version 15.1
 syntax [, n(integer 3000) sx(integer 0) s(integer 0) var(real 0.5) beta0(real -7.0) beta(real 0.5) kaC(real 1.25) ka(real 1.25)   ///
			a1(real 2.0) a0(real 2.0) gap(integer 30) c0(real -8.0) cx1(real 0.5) cx2(real 0.5) maxE(integer 100) ]

drop _all
local obs = `n'
set obs `obs'

gen id = _n

*follow-up time in days
gen tau = runiform(200,300)

*time constant covariates
gen z = rbinomial(1,0.5)
gen x = rbinomial(1,0.5)

if (`sx' != 0) {
	replace x = log(rgamma(4,1/4))
	}
	
if (`s' != 0) {
	replace tau = runiform(150,350)
	} 

gen exit = tau

*time to receive treatment
*hazard of treatment: Weibull disctribution with kaC>1
gen lamC = 1.0*exp(`c0'+`cx1'*x+0.5*z+rnormal(0,0.1))
*Underlying treatment time
gen ct0 = (-ln(runiform())/lamC)^(1/`kaC')

gen ct = ct0
replace ct = exit if ct>exit

*time to recurrent events
*hazard of events: Weibull disctribution with ka>1
local lbeta = log(`beta')

gen v = rnormal(0,sqrt(`var'))
gen lam = 0.5*exp(`beta0'+`cx2'*x+0.5*z+v)

gen t1 = (-ln(runiform())/lam)^(1/`ka')
replace t1 = (ct^`ka'-ln(runiform())/(exp(`lbeta'*(ct<exit))*lam))^(1/`ka')  if t1>ct 
replace t1 = exit if t1>exit

qui{
forvalues i = 2(1)`maxE'{
local j = `i'-1

gen ct_star = ((ct0^`kaC'-t`j'^`kaC')/`a1'+t`j'^`kaC')^(1/`kaC') if t`j'<exit & ct>t`j'
replace ct_star = ((ct0^`kaC'-(1-`a1')*t`j'^`kaC')/`a0'-(`a1'/`a0'-1)*(t`j'+`gap'/2)^`kaC' )^(1/`kaC')  if ct_star>t`j'+`gap'/2 & ct_star<.
replace ct_star = (ct0^`kaC'-(`a1'-`a0')*(t`j'+`gap'/2)^`kaC'-(`a0'-1)*(t`j'+`gap')^`kaC'+(`a1'-1)*t`j'^`kaC')^(1/`kaC')   if ct_star>t`j'+`gap' & ct_star<.

replace ct = ct_star  if ct_star<.
replace ct = exit  if ct>exit
replace ct0 = ct_star  if ct_star<. 

gen t`i' = (t`j'^`ka'-ln(runiform())/(exp(`lbeta'*(ct<exit)*(t`j'>ct))*lam))^(1/`ka')
replace t`i' = (ct^`ka'-ln(runiform())/(exp(`lbeta'*(ct<exit))*lam))^(1/`ka')  if t`i'>ct & t`j'<ct

replace t`i' = exit if t`i'>exit
drop ct_star
	}
}

drop ct0 lam* v tau
cap drop trt
gen trt = 0
replace trt = 1 if ct<exit

reshape long t ,i(id) j(temp)
drop temp 

duplicates drop
drop if t==0
 
gen event = 1
replace event = 0 if  t==exit 


by id (t),sort: gen n = _n
by id (t),sort: gen N = _N
tab event if n==N
gen ep = 1
replace ep = 2 if n==N
expand ep
drop n N
sort id t ep
by id t ep, sort: gen n = _n
by id t ep, sort: replace t = ct if n==2
drop n ep

drop if t>exit
duplicates drop

gen start = 0
gen end = t
by id (t), sort: replace start = end[_n-1] if _n>1
drop t

*time-varying variable 'treat'
gen treat=0
replace treat=1 if start>=ct

*observed time-constant covariate: 'x' 'z' 'z2'
********************************************
** AG model 1, Time-on-study: Calender time as timescale
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
by id: gen n = _n

reshape wide start end event, i(id) j(n)

gen cut = ct
replace cut = exit if cut>=exit

gen entry0 = 0
stset cut, failure(trt) enter(entry0) 

gen inclusive = (trt==0)
*match on z only:
qui sttocc_rev, match(z) n(1) nodots
gen index = ct
by _set (_case), sort: replace index = index[_N]

by _set, sort: gen Nset = _N
drop if Nset==1

drop Nset _set _time _case  cut inclusive entry0  

reshape long start end event, i(id) j(temp)
drop temp
drop if event==.

expand 2 if trt==0

sort id start end
by id start: gen n = _n
replace end = index if n==1 & end>index & start<index & trt==0
replace start = index if n==2 & end>index & start<index & trt==0
drop n
duplicates drop

replace event = 0 if end==index & trt==0

*data analysis from here:	
gen period = cond(end<=index,0,1)

cap drop n
by id (start),sort: gen n = _n

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


*PERR from AG model
gen pt = trt*period

stset end, failure(event) exit(t .) id(id) enter(start) 
tab _st,m
stcox  pt trt period , hr cluster(id)
matrix table_AG = r(table)
return scalar coef_AG = table_AG[1,1]
return scalar se_AG = table_AG[2,1]
return scalar pvalue_AG = table_AG[4,1]
return scalar cp_AG = cond(`beta'>=table_AG[5,1]&`beta'<=table_AG[6,1],1,0)

stset, clear

****************************
*find the estimator of `theta' and `delta' 
preserve

gen periodGap = cond(end<=index,0,2)

local m = 5
local mth = 10*`m'

qui forvalues k = 1(1)`m' {

expand 2
cap drop n
sort id start end
by id start: gen n = _n

local p = `mth'*`k'/`m'

replace event = 0 if  n==1 & end>index-`p' & start<index-`p'
replace end = index-`p' if n==1 & end>index-`p' & start<index-`p'
replace start = index-`p' if n==2 & end>index-`p' & start<index-`p'

drop n
duplicates drop

replace periodGap = 10+`k' if end>index-`mth'*`k'/`m' & end<=index-`mth'*(`k'-1)/`m'
}

*data analysis from here:	
gen ptGap = trt*periodGap
stset end, failure(event) exit(t .) id(id) enter(start) 

stcox  i.ptGap trt i.periodGap, nohr cluster(id)
matrix table_AG_Gap = r(table)

local gg = (table_AG_Gap[4,3]<0.05)
forvalues i = 2(1)4 {
local g = cond(`gg'==`i'-1,`gg'+(table_AG_Gap[4,2+`i']<0.05), `gg')
local gg = `g'
}

return scalar M = `g'

local p0 = `g'*10

stset, clear
restore

********************
foreach w in 1 2 {

preserve

local p = cond(`w'==1,`gap',`p0')

expand 2 

cap drop n
sort id start end
by id start: gen n = _n

replace event = 0 if  n==1 & end>index-`p' & start<index-`p'
replace end = index-`p' if n==1 & end>index-`p' & start<index-`p'
replace start = index-`p' if n==2 & end>index-`p' & start<index-`p'

drop n
duplicates drop

gen period2 = cond(end<=index,1,2)
replace period2 = 0 if end<=index-`p'

gen pt2 = trt*period2

stset end , failure(event) exit(t .) id(id) enter(start) 
stcox  i.pt2 trt i.period2,  cluster(id)

matrix table_AG`w' = r(table)
return scalar coef_AG`w' = table_AG`w'[1,3]
return scalar se_AG`w' = table_AG`w'[2,3]
return scalar pvalue_AG`w' = table_AG`w'[4,3]
return scalar cp_AG`w' = cond(`beta'>=table_AG`w'[5,3]&`beta'<=table_AG`w'[6,3],1,0)

return scalar coef_b3_AG`w' = table_AG`w'[1,2]
return scalar se_b3_AG`w' = table_AG`w'[2,2]

restore				 
}

local a2 = table_AG2[1,2]
drop pt* period*

***************************
*Adjustment from here:
*shifting according to all previous events before treatment time
foreach d in 1 2 {

preserve

cap drop n
sort id start end
by id: gen n = _n
su n
local max = r(max)-1

reshape wide start end event, i(id) j(n)
gen index0 = index

qui forvalues j = 1(1)`max'{
local i = `j'+1
if (`d'==2 ) {
gen ct_star = (index-end`j'+`a2'*end`j')/`a2'  if index>end`j' & end`j'<exit & trt==0 & event`j'==1
replace ct_star = index-`p0'*(`a2'-1)  if ct_star<. & ct_star>end`j'+`p0' & trt==0 & end`j'<. & event`j'==1
}
else {
gen ct_star = (index-end`j'+`a1'*end`j')/`a1'  if index>end`j' & end`j'<exit & trt==0 & event`j'==1
replace ct_star = (index-end`j'-`gap'/2*`a1')/`a0'+end`j'+`gap'/2  if ct_star<. & ct_star>end`j'+`gap'/2 & trt==0 & end`j'<. & event`j'==1

replace ct_star = index-`gap'*(`a1'/2+`a0'/2-1)  if ct_star<. & ct_star>end`j'+`gap' & trt==0 & end`j'<. & event`j'==1
}
replace index = ct_star if ct_star<.

replace index0 = ct_star if ct_star<.
replace index0 = exit if index0>exit
drop ct_star 

}
drop index

reshape long start end event, i(id) j(temp)
drop temp
drop if event==.

drop if index0>=exit

expand 2 if trt==0
cap drop n
sort id start end
by id start: gen n = _n
replace end = index0 if n==1 & end>index0 & start<index0 & trt==0
replace start = index0 if n==2 & end>index0 & start<index0 & trt==0
drop n
duplicates drop

replace event = 0 if end==index0 & trt==0
******
gen period = cond(end<=index0,0,1)
gen pt = trt*period

*PERR from AG model:
stset end, failure(event) exit(t .) id(id) enter(start) 
stcox  pt trt period,  cluster(id)

matrix table_AG_adj`d' = r(table)
return scalar coef_AG_adj_all`d' = table_AG_adj`d'[1,1]
return scalar se_AG_adj_all`d' = table_AG_adj`d'[2,1]
return scalar pvalue_AG_adj_all`d' = table_AG_adj`d'[4,1]
return scalar cp_AG_adj_all`d' = cond(`beta'>=table_AG_adj`d'[5,1]&`beta'<=table_AG_adj`d'[6,1],1,0)

restore				 
}

*shifting according to lastest event before treatment time

foreach d in 1 2 {

preserve

cap drop n
sort id start end
by id: gen n = _n
su n
local max = r(max)-1

reshape wide start end event, i(id) j(n)
gen index0 = index

qui forvalues j=1(1)`max'{
local i=`j'+1
if (`a0'>=1.0) {

if (`d'==2 ) {
gen ct_star = (index-end`j'+`a2'*end`j')/`a2'  if index>end`j' & index<=end`i'& end`j'<exit & trt==0 & event`j'==1
replace ct_star = index-`p0'*(`a2'-1)  if ct_star<. & ct_star>end`j'+`p0' & trt==0 & end`j'<. & event`j'==1
}
else {
gen ct_star = (index-end`j'+`a1'*end`j')/`a1'  if index>end`j' & index<=end`i'& end`j'<exit & trt==0 & event`j'==1
replace ct_star = (index-end`j'-`gap'/2*`a1')/`a0'+end`j'+`gap'/2  if ct_star<. & ct_star>end`j'+`gap'/2 & trt==0 & end`j'<. & event`j'==1

replace ct_star = index-`gap'*(`a1'/2+`a0'/2-1)  if ct_star<. & ct_star>end`j'+`gap' & trt==0 & end`j'<. & event`j'==1
}
	} 
else {

if (`d'==2 ) {
gen ct_star = (index-end`j'+`a2'*end`j')/`a2'  if index>end`j' & index<=end`i'& end`j'<exit & trt==0 & event`j'==1
replace ct_star = index-`p0'*(`a2'-1)  if ct_star<. & ct_star>end`j'+`p0' & trt==0 & end`j'<. & event`j'==1
replace ct_star = index+(end`i'-index)*0.99 if event`i'==1 & ct_star>end`i' & ct_star<.

}
else {
gen ct_star = (index-end`j'+`a1'*end`j')/`a1'  if index>end`j' & index<=end`i'& end`j'<exit & trt==0 & event`j'==1
replace ct_star = (index-end`j'-`gap'/2*`a1')/`a0'+end`j'+`gap'/2  if ct_star<. & ct_star>end`j'+`gap'/2 & trt==0 & end`j'<. & event`j'==1

replace ct_star = index-`gap'*(`a1'/2+`a0'/2-1)  if ct_star<. & ct_star>end`j'+`gap' & trt==0 & end`j'<. & event`j'==1
replace ct_star = index+(end`i'-index)*0.99 if event`i'==1 & ct_star>end`i' & ct_star<.

}
	}

replace index = ct_star if ct_star<.

replace index0 = ct_star if ct_star<.
replace index0 = exit if index0>exit
drop ct_star 

}
drop index

reshape long start end event, i(id) j(temp)
drop temp
drop if event==.

drop if index0>=exit

expand 2 if trt==0
cap drop n
sort id start end
by id start: gen n = _n
replace end = index0 if n==1 & end>index0 & start<index0 & trt==0
replace start = index0 if n==2 & end>index0 & start<index0 & trt==0
drop n
duplicates drop

replace event = 0 if end==index0 & trt==0

******
gen period = cond(end<=index0,0,1)
gen pt = trt*period

*PERR from AG model:
stset end , failure(event) exit(t .) id(id) enter(start) 
stcox  pt trt period ,  cluster(id)

matrix table_AG_adj`d' = r(table)
return scalar coef_AG_adj_last`d' = table_AG_adj`d'[1,1]
return scalar se_AG_adj_last`d' = table_AG_adj`d'[2,1]
return scalar pvalue_AG_adj_last`d' = table_AG_adj`d'[4,1]
return scalar cp_AG_adj_last`d' = cond(`beta'>=table_AG_adj`d'[5,1]&`beta'<=table_AG_adj`d'[6,1],1,0)

restore				 
}

eret clear

end

