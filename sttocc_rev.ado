*! version 6.0.4  21apr1999
** Convert survival-time data to case–control data
** original code from https://www.stata.com/updates/ado/sttocc.ado
** revise function sttocc in order to unique ID for recurrent event times
** Jul2024
** */ `include' & abs(`1' - `1'[`case'])<0.1    /* matching variable  */
** 	/* sorts with cases first, in order of timeout */


capture program drop sttocc_rev

program define sttocc_rev, rclass
	version 6
	st_is 2 analysis
	syntax [varlist] [, noDOTs GENerate(string) Match(varlist) /*
			*/  Number(integer 1) noSHow ]
	tokenize "`generate'"
	local d "`1'" 
	if "`d'"=="" {
		local d "_case"
	}
	local setid "`2'"
	if "`setid'"=="" {
		local setid "_set"
	}
	local ctime "`3'"
	if "`ctime'"=="" {
		local ctime "_time"
	}
	confirm new var `d' `setid' `ctime'
	if "`4'" != "" {
		di in red "more than 3 variables specified with generate option"
		di in red "only 3 variables are generated"
	        exit 198
	}
	st_show `show'
	if "`match'"!="" {
		di in gr "       matching for:  " in ye "`match'"
	}

	local wv : char _dta[st_wv]
	local ttype : type _t 
	local ftype : type _d
	
	* trim down to the dataset to be saved
	keep `varlist' _d  _t  _t0  `wv' `match'
tempvar mid
gen `mid'=uniform()

	tempvar idd sel include timein timeout new tied
	tempname nc obs tmpobs incr jitter

	* generate d, timein and timeout variables

	gen `ftype' `d' = _d 
	gen `ttype' `timeout' = _t 
	gen `ttype' `timein' = _t0 

	qui {

		/* checks on ties in timeout */
		/* place censorings after failures 
			and randomly jitter failures */
	
		summ `timeout'
		local jitter = sqrt(r(Var))*0.00005

		sort `timeout' `d' `mid'
		by `timeout' `d': gen int `tied' = (`d'!=0 & _N!=1)
		by `timeout' `d': replace `tied' = sum(`tied')
		by `timeout' `d': replace `tied' = 0 if _n<_N
		count if `tied'>0
		local nties = r(N)
		if `nties'>0 {
			noi di
			noi di in gr /*
			*/ "There were `nties' tied times involving failure(s)"
			noi di in gr _col(3) /*
			*/ "- tied failure times split at random"
			by `timeout': replace `timeout'=`timeout' + /*
			*/ cond(`d'==0, 0, uniform()/1000) /*
			*/ if _N != 1
		}  	
		/* revised by emily */

		/* sorts with cases first, in order of timeout */

		replace `d' = 1-`d'
		sort `d'  `mid'
		*sort `d' `timeout' `mid'
		replace `d' = 1-`d'
	
		/* counts number of cases */
	
		count if `d'==1
		scalar `nc'=r(N)
		noi di
		noi di in gr "There are " in ye `nc' in gr " cases"
	
		/* sets initial values */
	
		gen int `idd'=_n
		gen int `setid'=.
		scalar `obs'=_N
		scalar `tmpobs'=_N 
		gen byte `include'=.
		gen `ttype' `ctime'=.
	
		/* sample each set in turn */
	
		noi di in gr "Sampling " in ye "`number'" in gr /*
		*/ " controls for each case "
		local case = 1
		while `case'<=`nc' {
			if "`dots'"=="" {
				noi di  "." _continue
			}
			local ftime = `timeout'[`case']
		***	replace `include'= (_n != `case')
			replace `include'= (_n != `case' & inclusive==1 & id !=id[`case'] )      /* revised by emily */
			tokenize "`match'"
			while "`1'" != "" {
				/* sets include=1 for all records which pass
				the matching criteria except the current case */
				replace `include' = /*
				*/ `include' & abs(`1' - `1'[`case'])<0.1    /* matching variable  */
				mac shift
       		  	}
			replace `include'=`include' & `timein'<`ftime' /*
			*/ & `timeout'>`ftime' & `idd'!=.
	
			/* selects random sample from all records with
							 include = 1 */
			noi RSamp `idd' if `include', gen(`sel') n(`number')
			/* counts how many controls are selected */

			count if `sel'==1
			scalar `incr'=r(N) + 1
	
			/* expands selected controls and current case*/
	
			expand 2 if `idd'==`idd'[`case'] | `sel'==1
			replace `setid'=`case' if _n>`tmpobs'
			replace `ctime'=`timeout'[`case'] if _n>`tmpobs'
			replace `d'=0 if /*
			*/ `d'==1 & _n> `tmpobs' & `idd' != `idd'[`case']
			replace `idd'=. if _n>`obs'
			local case = `case'+1
			scalar `tmpobs' = `tmpobs' + `incr'
			* 13 Jan 2023: each person can only be selected once.
			replace inclusive=0 if `sel'==1
			drop `sel'
		}
		drop if _n<=`obs'
	}
	di
	sort `setid' `d' `mid'
	st_set clear
end

capture program drop RSamp

program define RSamp, rclass
	version 6 
	syntax varname [if] [, GENerate(string) Number(integer 1) ]
	tokenize `varlist'
	local id `1'
	tempvar u include
	confirm new var `generate'
	qui {
		qui gen `generate'=0
		qui count `if' `in'
		if r(N)<`number' {
			noi di in bl /*
			*/ " Warning: sample requested greater" /*
			*/ " than population. Only " r(N) " controls selected"
			replace `generate' = 1 `if' `in'
			exit 
		}
		qui gen `u'=uniform() `if' `in'
		sort `u'
		qui replace `generate'=1 in 1/`number'
		sort `id'
	}
end
