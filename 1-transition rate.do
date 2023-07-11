* Aim: Calculate the overall transition rates for subsequent estimation of disease-free LE

use main_99999.dta, clear
keep if _trans == 1
stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
stsplit age, every(1)
gen person_year = _t - _t0
preserve

****** empirical rates
foreach sex of numlist 0/1 {
    tab age _status if is_female == `sex', missing matrow(low_age) matcell(case)
    mat age_case = [low_age,case[1...,2..2]]
    loc interval_`sex' = rowsof(age_case)

    tabstat person_year if is_female == `sex', stat(sum) by(age) save
    mat null = r(Stat1)
    foreach j of numlist 2/`interval_`sex'' {
        mat null = [null\r(Stat`j')]
    }
    mat case_py_`sex' = [age_case,null]
}
mat case_py = [case_py_0\case_py_1]

clear
set obs 0
svmat case_py, names(var)
rename var1 age
rename var2 case1
rename var3 person_year1
gen sex = 1 if _n==1
    replace sex = 2 if _n !=1 & age[_n] < age[_n-1]
    replace sex = sex[_n-1] if age[_n] > age[_n-1]
gen empir_hazard1 = case1/person_year1
order sex age case1 person_year1 empir_hazard1

****** smoothed rates
/** The Gompertz model has been fitted with the Poisson distribution for the number of deaths and a log link to ages.
ref: https://www.longevitas.co.uk/site/informationmatrix/anotherlookatthegompertzmodel.html**/
restore
sort is_female
gen nage_M = ${b_age} if is_female == 0 & _n == 1
    replace nage_M = nage_M[_n-1] + 1 if is_female == 0 & _n > 1
gsort -is_female
gen nage_F = ${b_age} if is_female == 1 & _n == 1
    replace nage_F = nage_F[_n-1] + 1 if is_female == 1 & _n > 1
egen nage = rowmax(nage_M nage_F)
drop nage_M nage_F

* male
poisson _d age if is_female == 0, exposure(person_year)
mat A = r(table)
loc b = A[1,1]
loc a = A[1,2]
gen hazard1_M = exp(`b'*nage + `a') if is_female == 0
* female
poisson _d age if is_female == 1, exposure(person_year)
mat A = r(table)
loc b = A[1,1]
loc a = A[1,2]
gen hazard1_F = exp(`b'*nage + `a') if is_female == 1

egen hazard1 = rowmax(hazard1_M hazard1_F)
keep is_female nage hazard1
rename is_female sex
    replace sex = sex + 1
drop if nage < ${b_age} | nage > ${e_age}
sort sex nage
rename nage age
