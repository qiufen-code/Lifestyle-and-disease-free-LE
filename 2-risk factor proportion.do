*** Aim: sex- and age-specific risk factor prevalence for each transition

cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
use main_99999.dta, clear
keep if _trans == 1
stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
stsplit age, every(5)
    replace age=80 if age>80
    
* male
sum age if is_female==0
loc min = r(min)
loc max = r(max)
levelsof HLS_5groups if is_female==0, loc(level) clean

proportion HLS_5groups if is_female==0, over(age) cformat(%9.3f)
mat A=r(table)
foreach age of numlist `min'(5)`max' {
    foreach l of loc level {
        loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
        loc p`age'_`l'=A[1,`col']
        loc se`age'_`l'=A[2,`col']
    }
}

* female
sum age if is_female==1
loc min = r(min)
loc max = r(max)
levelsof HLS_5groups if is_female==1, loc(level) clean

proportion HLS_5groups if is_female==1, over(age) cformat(%9.3f)
mat A=r(table)
foreach age of numlist `min'(5)`max' {
    foreach l of loc level {
        loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
        loc p`age'_`l'=A[1,`col']
        loc se`age'_`l'=A[2,`col']
    }
}
