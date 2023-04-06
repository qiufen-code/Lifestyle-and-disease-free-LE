* Aim: Calculate the overall transition rates for subsequent estimation of disease-free LE

*** define the start and end age
global b_age = 30
global e_age = 90

*** specify the database for use
loc prefix main

*** main code
foreach datas of loc prefix {

	* 在主分析中关注单一疾病及其组合，在敏感性分析或是亚组分析中只关注疾病组合
	if "`datas'" == "main" {
		loc incidence 99999 16418 16419 16420
	}
	if "`datas'" == "T2D" {
		loc incidence 00000
	}
	else {
		loc incidence 99999
	}


	log using "D:\SunQF\03-healthy_life_expectancy\Check\transition rate_`datas'.log", replace

	postfile se_rate sex disease trans coef cons se_coef se_cons using "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_rate_`datas'.dta", replace

	foreach i of loc incidence {
		display _dup(30) "*" "`datas'_`i'" _dup(30) "*"
		******************************* 1 disease incidence *******************************
		display _dup(30) "*" "1 disease incidence" _dup(30) "*"
		cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
		use `datas'_`i'.dta, clear
		keep if _trans == 1
		stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
		stsplit age, every(1)    //对数据集进行划分，得到attained age
		gen person_year = _t - _t0
		preserve
		
		****** empirical rates
		display "* empirical rates"
		foreach sex of numlist 0/1 {
			display "* is_female = `sex'"
			tab age _status if is_female == `sex', missing matrow(low_age) matcell(case)
			mat age_case = [low_age,case[1...,2..2]]
			loc interval_`sex' = rowsof(age_case)      //计数分析人群随访期间的年龄跨度，例如30-91岁

			tabstat person_year if is_female == `sex', stat(sum) by(age) save
			mat null = r(Stat1)
			foreach j of numlist 2/`interval_`sex'' {
				mat null = [null\r(Stat`j')]
			}
			mat case_py_`sex' = [age_case,null]       //得到case和人时数据集
		}
		mat case_py = [case_py_0\case_py_1]           //将男性和女性的数据集合并，变量分别为age、case、PYs
		
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
		tempfile observed_rate`i'_1
		save "`observed_rate`i'_1'", replace
		
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
		display "* smoothed rates: male"
		poisson _d age if is_female == 0, exposure(person_year)
		mat A = r(table)
		post se_rate (1) (`i') (1) (A[1,1]) (A[1,2]) (A[2,1]) (A[2,2])
		
		loc b = A[1,1]
		loc a = A[1,2]
		gen hazard1_M = exp(`b'*nage + `a') if is_female == 0
		* female
		display "* smoothed rates: female"
		poisson _d age if is_female == 1, exposure(person_year)
		mat A = r(table)
		post se_rate (2) (`i') (1) (A[1,1]) (A[1,2]) (A[2,1]) (A[2,2])
		
		loc b = A[1,1]
		loc a = A[1,2]
		gen hazard1_F = exp(`b'*nage + `a') if is_female == 1

		egen hazard1 = rowmax(hazard1_M hazard1_F)
		keep is_female nage hazard1
		rename is_female sex
			replace sex = sex + 1                 //更改is_female(0=men, 1=women)变量为sex(1=men, 2=women)
		drop if nage < ${b_age} | nage > ${e_age}
		sort sex nage
		rename nage age
		
		tempfile predict_rate`i'_1
		save "`predict_rate`i'_1'", replace
		
		******************************* 2 die without disease *******************************
		display _dup(30) "*" "2 die without disease" _dup(30) "*"
		cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
		use `datas'_`i'.dta, clear
		keep if _trans == 2
		stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
		stsplit age, every(1)
		gen person_year = _t - _t0
		preserve
		
		****** empirical rates
		display "* empirical rates"
		foreach sex of numlist 0/1 {
			display "* is_female = `sex'"
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
		rename var2 case2
		rename var3 person_year2
		gen sex = 1 if _n==1
			replace sex = 2 if _n !=1 & age[_n] < age[_n-1]
			replace sex = sex[_n-1] if age[_n] > age[_n-1]
		gen empir_hazard2 = case2/person_year2
		order sex age case2 person_year2 empir_hazard2
		tempfile observed_rate`i'_2
		save "`observed_rate`i'_2'", replace
		
		****** smoothed rates
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
		display "* smoothed rates: male"
		poisson _d age if is_female == 0, exposure(person_year)
		mat A = r(table)
		post se_rate (1) (`i') (2) (A[1,1]) (A[1,2]) (A[2,1]) (A[2,2])
		
		loc b = A[1,1]
		loc a = A[1,2]
		gen hazard2_M = exp(`b'*nage + `a') if is_female == 0
		* female
		display "* smoothed rates: female"
		poisson _d age if is_female == 1, exposure(person_year)
		mat A = r(table)
		post se_rate (2) (`i') (2) (A[1,1]) (A[1,2]) (A[2,1]) (A[2,2])
		
		loc b = A[1,1]
		loc a = A[1,2]
		gen hazard2_F = exp(`b'*nage + `a') if is_female == 1

		egen hazard2 = rowmax(hazard2_M hazard2_F)
		keep is_female nage hazard2
		rename is_female sex
			replace sex = sex + 1
		drop if nage < ${b_age} | nage > ${e_age}
		sort sex nage
		rename nage age
		
		tempfile predict_rate`i'_2
		save "`predict_rate`i'_2'", replace

		******************************* 3 die after disease *******************************
		dis _dup(30) "*" "3 die after disease" _dup(30) "*"
		cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
		use `datas'_`i'.dta, clear
		keep if _trans == 3
		stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
		stsplit age, every(1)
		* list ep_`i'_combined_ep ep_16416_du_ep ep_`i'_combined_ep_date ep_16416_du_ep_date censoring_reason if age==.
			clonevar age_miss = age
				replace age=floor((_start-dob_anon)/365.25) if age_miss==.
				replace _t0=((_start-0.5)-dob_anon)/365.25 if age_miss==.     //对于在随访结束的当天发病或发病当天立即失访者，将发病时间减去0.5天
				replace _t=(_stop-dob_anon)/365.25 if age_miss==.
				replace _d=0 if age_miss==.                                   //发生上述两种情形者，都没有发生病后死亡
		gen person_year = _t - _t0
		preserve
		
		****** empirical rates
		display "* empirical rates"
		foreach sex of numlist 0/1 {
			display "* is_female = `sex'"
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
		rename var2 case3
		rename var3 person_year3
		gen sex = 1 if _n==1
			replace sex = 2 if _n !=1 & age[_n] < age[_n-1]
			replace sex = sex[_n-1] if age[_n] > age[_n-1]
		gen empir_hazard3 = case3/person_year3
		order sex age case3 person_year3 empir_hazard3
		tempfile observed_rate`i'_3
		save "`observed_rate`i'_3'", replace
		
		****** smoothed rates
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
		display "* smoothed rates: male"
		poisson _d age if is_female == 0, exposure(person_year)
		mat A = r(table)
		post se_rate (1) (`i') (3) (A[1,1]) (A[1,2]) (A[2,1]) (A[2,2])
		
		loc b = A[1,1]
		loc a = A[1,2]
		gen hazard3_M = exp(`b'*nage + `a') if is_female == 0
		* female
		display "* smoothed rates: female"
		poisson _d age if is_female == 1, exposure(person_year)
		mat A = r(table)
		post se_rate (2) (`i') (3) (A[1,1]) (A[1,2]) (A[2,1]) (A[2,2])
		
		loc b = A[1,1]
		loc a = A[1,2]
		gen hazard3_F = exp(`b'*nage + `a') if is_female == 1

		egen hazard3 = rowmax(hazard3_M hazard3_F)
		keep is_female nage hazard3
		rename is_female sex
			replace sex = sex + 1
		drop if nage < ${b_age} | nage > ${e_age}
		sort sex nage
		rename nage age

		tempfile predict_rate`i'_3
		save "`predict_rate`i'_3'", replace

		* merge observed and predicted rate
		use "`observed_rate`i'_1'", clear
		merge 1:1 sex age using "`observed_rate`i'_2'", nogenerate
		merge 1:1 sex age using "`observed_rate`i'_3'", nogenerate
		merge 1:1 sex age using "`predict_rate`i'_1'", nogenerate
		merge 1:1 sex age using "`predict_rate`i'_2'", nogenerate
		merge 1:1 sex age using "`predict_rate`i'_3'", nogenerate
		order sex age case1 person_year1 empir_hazard1 case2 person_year2 empir_hazard2 case3 person_year3 empir_hazard3 hazard1 hazard2 hazard3
		save D:\SunQF\03-healthy_life_expectancy\Check\trans_rate\rate`i'_`datas'.dta, replace
		export excel using "D:\SunQF\03-healthy_life_expectancy\Check\trans_rate\rate`i'_`datas'.xlsx", firstrow(variables) replace
	}
	postclose se_rate
	log close

	* For CI
	use "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_rate_`datas'.dta", clear
	export excel using "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_rate_`datas'.xlsx", firstrow(variables) replace
	erase "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_rate_`datas'.dta"
}

