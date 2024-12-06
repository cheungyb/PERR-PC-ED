
capture program drop Cox_Ratio

*Cox model to first event

program define Cox_Ratio, rclass
version 15.1

args newpairid trt start end period  post first1  prior first0

*use only data from the post-period 
preserve
keep if `period'==1
stset `end' if `post'<=`first1', failure(`first1')   id(`newpairid') enter(`start') 
stcox `trt', nohr
scalar Cox1_b1=_b[`trt']
restore

*use only data from the prior-period
preserve
keep if `period'==0
stset `end' if `prior'<=`first0' , failure(`first0')  id(`newpairid') enter(`start') 
stcox `trt', nohr
scalar Cox0_b1=_b[`trt']
restore

return scalar ratio = Cox1_b1-Cox0_b1

end

