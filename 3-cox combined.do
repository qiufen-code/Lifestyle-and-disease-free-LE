* Cox proportional hazards analysis

* (1) Incidence
cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
use main_99999.dta, clear
keep if _trans == 1
stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
gen person_year = _t - _t0

* male
count if is_female == 0 & _status == 1
loc case = r(N)
tabstat person_year if is_female == 0, stat(sum) save
mat a=r(StatTotal)
loc py = a[1,1]
loc rate = `case'/`py'*1000
loc rate_1dec = string(`rate',"%9.1f")

* Person years & events
tab HLS_5groups if is_female == 0 & _status == 1
tabstat person_year if is_female == 0, stat(sum) by(HLS_5groups)

display "* Model1: null model"
stcox i.HLS_5groups if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
display "* Model2: Model1 + sociodemographic factors + family history"
stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 

* female
count if is_female == 1 & _status == 1
loc case = r(N)
tabstat person_year if is_female == 1, stat(sum) save
mat a=r(StatTotal)
loc py = a[1,1]
loc rate = `case'/`py'*1000
loc rate_1dec = string(`rate',"%9.1f")

* Person years & events
tab HLS_5groups if is_female == 1 & _status == 1
tabstat person_year if is_female == 1, stat(sum) by(HLS_5groups)

display "* Model1: null model"
stcox i.HLS_5groups if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
display "* Model2: Model1 + sociodemographic factors + family history + reproductive history"
stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 

* (2) death without disease 
cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
use main_99999.dta, clear
keep if _trans == 2
stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
gen person_year = _t - _t0

* male
count if is_female == 0 & _status == 1
loc case = r(N)
tabstat person_year if is_female == 0, stat(sum) save
mat a=r(StatTotal)
loc py = a[1,1]
loc rate = `case'/`py'*1000
loc rate_1dec = string(`rate',"%9.1f")

* Person years & events
tab HLS_5groups if is_female == 0 & _status == 1
tabstat person_year if is_female == 0, stat(sum) by(HLS_5groups)

display "* Model1: null model"
stcox i.HLS_5groups if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
display "* Model2: Model1 + sociodemographic factors + family history"
stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 

* female
count if is_female == 1 & _status == 1
loc case = r(N)
tabstat person_year if is_female == 1, stat(sum) save
mat a=r(StatTotal)
loc py = a[1,1]
loc rate = `case'/`py'*1000
loc rate_1dec = string(`rate',"%9.1f")

* Person years & events
tab HLS_5groups if is_female == 1 & _status == 1
tabstat person_year if is_female == 1, stat(sum) by(HLS_5groups)

display "* Model1: null model"
stcox i.HLS_5groups if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
display "* Model2: Model1 + sociodemographic factors + family history + reproductive history"
stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 

* (3) death after disease
cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
use main_99999.dta, clear
keep if _trans == 3
replace _start = _start - 0.5 if _start == _stop
stset _stop, id(csid) enter(_start) origin(time _start) failure(_status==1) scale(365.25)
gen person_year = _t - _t0

* male
count if is_female == 0 & _status == 1
loc case = r(N)
tabstat person_year if is_female == 0, stat(sum) save
mat a=r(StatTotal)
loc py = a[1,1]
loc rate = `case'/`py'*1000
loc rate_1dec = string(`rate',"%9.1f")

* Person years & events
tab HLS_5groups if is_female == 0 & _status == 1
tabstat person_year if is_female == 0, stat(sum) by(HLS_5groups)

display "* Model1: null model "
stcox i.HLS_5groups if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
display "* Model2: Model1 + diagnosis age"
stcox i.HLS_5groups age_at_99999_diagnosis if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
display "* Model3: Model2 + sociodemographic factors + family history"
stcox i.HLS_5groups age_at_99999_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 

* female
count if is_female == 1 & _status == 1
loc case = r(N)
tabstat person_year if is_female == 1, stat(sum) save
mat a=r(StatTotal)
loc py = a[1,1]
loc rate = `case'/`py'*1000
loc rate_1dec = string(`rate',"%9.1f")

* Person years & events
tab HLS_5groups if is_female == 1 & _status == 1
tabstat person_year if is_female == 1, stat(sum) by(HLS_5groups)

display "* Model1: null model "
stcox i.HLS_5groups if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
display "* Model2: Model1 + diagnosis age"
stcox i.HLS_5groups age_at_99999_diagnosis if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
display "* Model3: Model2 + sociodemographic factors + family history + reproductive history"
stcox i.HLS_5groups age_at_99999_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
