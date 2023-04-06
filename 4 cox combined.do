* Cox proportional hazards analysis

*** define the family history variable for each transition
global family_history_99999 "i.heart_attack_FH i.stroke_FH i.cancer_FH"                  //CVD、癌症、CRD
global family_history_16418 "i.heart_attack_FH i.stroke_FH"
global family_history_16419 "i.cancer_FH"
global family_history_16420
global family_history_16421 "i.diabetes_FH"                                              //T2D
global family_history_00000 "i.heart_attack_FH i.stroke_FH i.cancer_FH i.diabetes_FH"    //CVD、癌症、CRD和T2D

*** define the start and end age
global b_age = 30
global e_age = 90
loc prefix main

*** main code
foreach datas of loc prefix {

	* 在主分析中关注单一疾病及其组合，在敏感性分析或是亚组分析中只关注疾病组合
	if "`datas'" == "main" | "`datas'" == "exclu2y" {
		loc incidence 99999 16418 16419 16420
	}
	else if "`datas'" == "T2D" {
		loc incidence 00000
	}
	else {
		loc incidence 99999
	}


	/* combined lifestyle */
	log using "D:\SunQF\03-healthy_life_expectancy\Check\cox_combined_`datas'.log", replace

	* (1) Incidence
	display _dup(30) "*" "incidence" _dup(30) "*"

	postfile table1 str20 (sex outcome trans) cases rates str20 (hr1 hr2 hr3 hr4 hr5 trend pvalue) using "D:\SunQF\03-healthy_life_expectancy\3 分析结果\Table-1_`datas'.dta", replace
	postfile combined_trans1 sex age incidence hr11 hr12 hr13 hr14 hr15 using "D:\SunQF\03-healthy_life_expectancy\Check\Cox\combined_trans1.dta", replace
	postfile se str20 (sex life) incidence trans beta1 beta2 beta3 beta4 beta5 se1 se2 se3 se4 se5 using "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_cox_combined_`datas'.dta", replace

	foreach i of loc incidence {
		display _dup(30) "*" "`datas'_`i'" _dup(30) "*"
		cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
		use `datas'_`i'.dta, clear
		keep if _trans == 1
		stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
		gen person_year = _t - _t0

		* male
		display _dup(30) "*" "Male" _dup(30) "*"
		display "* Person year & cases"
		count if is_female == 0 & _status == 1
		loc case = r(N)
		tabstat person_year if is_female == 0, stat(sum) save
		mat a=r(StatTotal)
		loc py = a[1,1]
		loc rate = `case'/`py'*1000               //每1000 PY的transition rate
		loc rate_1dec = string(`rate',"%9.1f")    //对transition rate保留1位小数

		* 查看每一个暴露亚组内的发生事件数，用于查证CI尤其宽的结果
		tab HLS_5groups if is_female == 0 & _status == 1
		tabstat person_year if is_female == 0, stat(sum) by(HLS_5groups)

		display "* Model1: null model"
		stcox i.HLS_5groups if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		display "* Model2: Model1 + sociodemographic factors + family history"
		stcox i.HLS_5groups i.highest_education i.marital_status ${family_history_`i'} if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		matrix A=r(table)
		foreach j of numlist ${b_age}/${e_age} {
			post combined_trans1 (1) (`j') (`i') (A[1,1]) (A[1,2]) (A[1,3]) (A[1,4]) (A[1,5])
		}

		* HR的自然对数beta及其标准误
		stcox i.HLS_5groups i.highest_education i.marital_status ${family_history_`i'} if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) nohr
		matrix u=r(table)
		post se ("male") ("HLS") (`i') (1) (u[1,1]) (u[1,2]) (u[1,3]) (u[1,4]) (u[1,5]) (u[2,1]) (u[2,2]) (u[2,3]) (u[2,4]) (u[2,5])

		* Table 1.
		loc cat_1=string(A[1,1],"%9.2f")+" ("+"Referent"+")"
		loc cat_2=string(A[1,2],"%9.2f")+" ("+string(A[5,2],"%9.2f")+"-"+string(A[6,2],"%9.2f")+")"
		loc cat_3=string(A[1,3],"%9.2f")+" ("+string(A[5,3],"%9.2f")+"-"+string(A[6,3],"%9.2f")+")"
		loc cat_4=string(A[1,4],"%9.2f")+" ("+string(A[5,4],"%9.2f")+"-"+string(A[6,4],"%9.2f")+")"
		loc cat_5=string(A[1,5],"%9.2f")+" ("+string(A[5,5],"%9.2f")+"-"+string(A[6,5],"%9.2f")+")"

		display "* Trend test"
		stcox c.HLS_5groups i.highest_education i.marital_status ${family_history_`i'} if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		matrix B=r(table)
		loc trend=string(B[1,1],"%9.2f")+" ("+string(B[5,1],"%9.2f")+"-"+string(B[6,1],"%9.2f")+")"
		loc p=string(B[4,1],"%9.2e")
		post table1 ("Male") ("`i'") ("Baseline → disease") (`case') (`rate_1dec') ("`cat_1'") ("`cat_2'") ("`cat_3'") ("`cat_4'") ("`cat_5'") ("`trend'") ("`p'")
		
		* female
		display _dup(30) "*" "Female" _dup(30) "*"
		display "* Person year & cases"
		count if is_female == 1 & _status == 1
		loc case = r(N)
		tabstat person_year if is_female == 1, stat(sum) save
		mat a=r(StatTotal)
		loc py = a[1,1]
		loc rate = `case'/`py'*1000
		loc rate_1dec = string(`rate',"%9.1f")

		* 查看每一个暴露亚组内的发生事件数
		tab HLS_5groups if is_female == 1 & _status == 1
		tabstat person_year if is_female == 1, stat(sum) by(HLS_5groups)

		display "* Model1: null model"
		stcox i.HLS_5groups if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		display "* Model2: Model1 + sociodemographic factors + family history + reproductive history"
		stcox i.HLS_5groups i.highest_education i.marital_status ${family_history_`i'} i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		matrix A=r(table)
		foreach j of numlist ${b_age}/${e_age} {
			post combined_trans1 (2) (`j') (`i') (A[1,1]) (A[1,2]) (A[1,3]) (A[1,4]) (A[1,5])
		}

		* HR的自然对数beta及其标准误
		stcox i.HLS_5groups i.highest_education i.marital_status ${family_history_`i'} i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) nohr
		matrix u=r(table)
		post se ("female") ("HLS") (`i') (1) (u[1,1]) (u[1,2]) (u[1,3]) (u[1,4]) (u[1,5]) (u[2,1]) (u[2,2]) (u[2,3]) (u[2,4]) (u[2,5])

		* Table 1.
		loc cat_1=string(A[1,1],"%9.2f")+" ("+"Referent"+")"
		loc cat_2=string(A[1,2],"%9.2f")+" ("+string(A[5,2],"%9.2f")+"-"+string(A[6,2],"%9.2f")+")"
		loc cat_3=string(A[1,3],"%9.2f")+" ("+string(A[5,3],"%9.2f")+"-"+string(A[6,3],"%9.2f")+")"
		loc cat_4=string(A[1,4],"%9.2f")+" ("+string(A[5,4],"%9.2f")+"-"+string(A[6,4],"%9.2f")+")"
		loc cat_5=string(A[1,5],"%9.2f")+" ("+string(A[5,5],"%9.2f")+"-"+string(A[6,5],"%9.2f")+")"
		
		display "* Trend test"
		stcox c.HLS_5groups i.highest_education i.marital_status ${family_history_`i'} i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		matrix B=r(table)
		loc trend=string(B[1,1],"%9.2f")+" ("+string(B[5,1],"%9.2f")+"-"+string(B[6,1],"%9.2f")+")"
		loc p=string(B[4,1],"%9.2e")
		post table1 ("Female") ("`i'") ("Baseline → disease") (`case') (`rate_1dec') ("`cat_1'") ("`cat_2'") ("`cat_3'") ("`cat_4'") ("`cat_5'") ("`trend'") ("`p'")
	}
	postclose combined_trans1

	* (2) death without disease 
	display _dup(30) "*" "death without disease" _dup(30) "*"

	postfile combined_trans2 sex age incidence hr21 hr22 hr23 hr24 hr25 using "D:\SunQF\03-healthy_life_expectancy\Check\Cox\combined_trans2.dta", replace

	foreach i of loc incidence {
		display _dup(30) "*" "`datas'_`i'" _dup(30) "*"
		cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
		use `datas'_`i'.dta, clear
		keep if _trans == 2
		stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
		gen person_year = _t - _t0

		* male
		display _dup(30) "*" "Male" _dup(30) "*"
		display "* Person year & cases"
		count if is_female == 0 & _status == 1
		loc case = r(N)
		tabstat person_year if is_female == 0, stat(sum) save
		mat a=r(StatTotal)
		loc py = a[1,1]
		loc rate = `case'/`py'*1000
		loc rate_1dec = string(`rate',"%9.1f")

		* 查看每一个暴露亚组内的发生事件数
		tab HLS_5groups if is_female == 0 & _status == 1
		tabstat person_year if is_female == 0, stat(sum) by(HLS_5groups)

		display "* Model1: null model"
		stcox i.HLS_5groups if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		display "* Model2: Model1 + sociodemographic factors + family history"
		stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		matrix A=r(table)
		foreach j of numlist ${b_age}/${e_age} {
			post combined_trans2 (1) (`j') (`i') (A[1,1]) (A[1,2]) (A[1,3]) (A[1,4]) (A[1,5])
		}

		* HR的自然对数beta及其标准误
		stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) nohr
		matrix u=r(table)
		post se ("male") ("HLS") (`i') (2) (u[1,1]) (u[1,2]) (u[1,3]) (u[1,4]) (u[1,5]) (u[2,1]) (u[2,2]) (u[2,3]) (u[2,4]) (u[2,5])

		* Table 1.
		loc cat_1=string(A[1,1],"%9.2f")+" ("+"Referent"+")"
		loc cat_2=string(A[1,2],"%9.2f")+" ("+string(A[5,2],"%9.2f")+"-"+string(A[6,2],"%9.2f")+")"
		loc cat_3=string(A[1,3],"%9.2f")+" ("+string(A[5,3],"%9.2f")+"-"+string(A[6,3],"%9.2f")+")"
		loc cat_4=string(A[1,4],"%9.2f")+" ("+string(A[5,4],"%9.2f")+"-"+string(A[6,4],"%9.2f")+")"
		loc cat_5=string(A[1,5],"%9.2f")+" ("+string(A[5,5],"%9.2f")+"-"+string(A[6,5],"%9.2f")+")"
		
		display "* Trend test"
		stcox c.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		matrix B=r(table)
		loc trend=string(B[1,1],"%9.2f")+" ("+string(B[5,1],"%9.2f")+"-"+string(B[6,1],"%9.2f")+")"
		loc p=string(B[4,1],"%9.2e")
		post table1 ("Male") ("`i'") ("Baseline → death") (`case') (`rate_1dec') ("`cat_1'") ("`cat_2'") ("`cat_3'") ("`cat_4'") ("`cat_5'") ("`trend'") ("`p'")
		
		* female
		display _dup(30) "*" "Female" _dup(30) "*"
		display "* Person year & cases"
		count if is_female == 1 & _status == 1
		loc case = r(N)
		tabstat person_year if is_female == 1, stat(sum) save
		mat a=r(StatTotal)
		loc py = a[1,1]
		loc rate = `case'/`py'*1000
		loc rate_1dec = string(`rate',"%9.1f")

		* 查看每一个暴露亚组内的发生事件数
		tab HLS_5groups if is_female == 1 & _status == 1
		tabstat person_year if is_female == 1, stat(sum) by(HLS_5groups)

		display "* Model1: null model"
		stcox i.HLS_5groups if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		display "* Model2: Model1 + sociodemographic factors + family history + reproductive history"
		stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		matrix A=r(table)
		foreach j of numlist ${b_age}/${e_age} {
			post combined_trans2 (2) (`j') (`i') (A[1,1]) (A[1,2]) (A[1,3]) (A[1,4]) (A[1,5])
		}

		* HR的自然对数beta及其标准误
		stcox i.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) nohr
		matrix u=r(table)
		post se ("female") ("HLS") (`i') (2) (u[1,1]) (u[1,2]) (u[1,3]) (u[1,4]) (u[1,5]) (u[2,1]) (u[2,2]) (u[2,3]) (u[2,4]) (u[2,5])

		* Table 1.
		loc cat_1=string(A[1,1],"%9.2f")+" ("+"Referent"+")"
		loc cat_2=string(A[1,2],"%9.2f")+" ("+string(A[5,2],"%9.2f")+"-"+string(A[6,2],"%9.2f")+")"
		loc cat_3=string(A[1,3],"%9.2f")+" ("+string(A[5,3],"%9.2f")+"-"+string(A[6,3],"%9.2f")+")"
		loc cat_4=string(A[1,4],"%9.2f")+" ("+string(A[5,4],"%9.2f")+"-"+string(A[6,4],"%9.2f")+")"
		loc cat_5=string(A[1,5],"%9.2f")+" ("+string(A[5,5],"%9.2f")+"-"+string(A[6,5],"%9.2f")+")"

		display "* Trend test"
		stcox c.HLS_5groups i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		matrix B=r(table)
		loc trend=string(B[1,1],"%9.2f")+" ("+string(B[5,1],"%9.2f")+"-"+string(B[6,1],"%9.2f")+")"
		loc p=string(B[4,1],"%9.2e")
		post table1 ("Female") ("`i'") ("Baseline → death") (`case') (`rate_1dec') ("`cat_1'") ("`cat_2'") ("`cat_3'") ("`cat_4'") ("`cat_5'") ("`trend'") ("`p'")
	}
	postclose combined_trans2

	* (3) death after disease
	display _dup(30) "*" "death without disease" _dup(30) "*"

	postfile combined_trans3 sex age incidence hr31 hr32 hr33 hr34 hr35 using "D:\SunQF\03-healthy_life_expectancy\Check\Cox\combined_trans3.dta", replace

	foreach i of loc incidence {
		display _dup(30) "*" "`datas'_`i'" _dup(30) "*"
		cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
		use `datas'_`i'.dta, clear
		keep if _trans == 3
		replace _start = _start - 0.5 if _start == _stop                                             //发病时间与失访或是随访终点时间重合者，将其发病时间减去0.5天
		stset _stop, id(csid) enter(_start) origin(time _start) failure(_status==1) scale(365.25)    //extended semi-markov模型
		gen person_year = _t - _t0

		* male
		display _dup(30) "*" "Male" _dup(30) "*"
		display "* Person year & cases"
		count if is_female == 0 & _status == 1
		loc case = r(N)
		tabstat person_year if is_female == 0, stat(sum) save
		mat a=r(StatTotal)
		loc py = a[1,1]
		loc rate = `case'/`py'*1000
		loc rate_1dec = string(`rate',"%9.1f")

		* 查看每一个暴露亚组内的发生事件数
		tab HLS_5groups if is_female == 0 & _status == 1
		tabstat person_year if is_female == 0, stat(sum) by(HLS_5groups)

		display "* Model1: null model "
		stcox i.HLS_5groups if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		display "* Model2: Model1 + 诊断年龄"
		stcox i.HLS_5groups age_at_`i'_diagnosis if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		display "* Model3: Model2 + sociodemographic factors + family history"
		stcox i.HLS_5groups age_at_`i'_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		matrix A=r(table)
		foreach j of numlist ${b_age}/${e_age} {
			post combined_trans3 (1) (`j') (`i') (A[1,1]) (A[1,2]) (A[1,3]) (A[1,4]) (A[1,5])
		}

		* HR的自然对数beta及其标准误
		stcox i.HLS_5groups age_at_`i'_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) nohr
		matrix u=r(table)
		post se ("male") ("HLS") (`i') (3) (u[1,1]) (u[1,2]) (u[1,3]) (u[1,4]) (u[1,5]) (u[2,1]) (u[2,2]) (u[2,3]) (u[2,4]) (u[2,5])

		* Table 1.
		loc cat_1=string(A[1,1],"%9.2f")+" ("+"Referent"+")"
		loc cat_2=string(A[1,2],"%9.2f")+" ("+string(A[5,2],"%9.2f")+"-"+string(A[6,2],"%9.2f")+")"
		loc cat_3=string(A[1,3],"%9.2f")+" ("+string(A[5,3],"%9.2f")+"-"+string(A[6,3],"%9.2f")+")"
		loc cat_4=string(A[1,4],"%9.2f")+" ("+string(A[5,4],"%9.2f")+"-"+string(A[6,4],"%9.2f")+")"
		loc cat_5=string(A[1,5],"%9.2f")+" ("+string(A[5,5],"%9.2f")+"-"+string(A[6,5],"%9.2f")+")"
		
		display "* Trend test"
		stcox c.HLS_5groups age_at_`i'_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH if is_female == 0, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		matrix B=r(table)
		loc trend=string(B[1,1],"%9.2f")+" ("+string(B[5,1],"%9.2f")+"-"+string(B[6,1],"%9.2f")+")"
		loc p=string(B[4,1],"%9.2e")
		post table1 ("Male") ("`i'") ("Disease → death") (`case') (`rate_1dec') ("`cat_1'") ("`cat_2'") ("`cat_3'") ("`cat_4'") ("`cat_5'") ("`trend'") ("`p'")
		
		* female
		display _dup(30) "*" "Female" _dup(30) "*"
		display "* Person year & cases"
		count if is_female == 1 & _status == 1
		loc case = r(N)
		tabstat person_year if is_female == 1, stat(sum) save
		mat a=r(StatTotal)
		loc py = a[1,1]
		loc rate = `case'/`py'*1000
		loc rate_1dec = string(`rate',"%9.1f")

		* 查看每一个暴露亚组内的发生事件数
		tab HLS_5groups if is_female == 1 & _status == 1
		tabstat person_year if is_female == 1, stat(sum) by(HLS_5groups)

		display "* Model1: null model "
		stcox i.HLS_5groups if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		display "* Model2: Model1 + 诊断年龄"
		stcox i.HLS_5groups age_at_`i'_diagnosis if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		display "* Model3: Model2 + sociodemographic factors + family history + reproductive history"
		stcox i.HLS_5groups age_at_`i'_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) 
		matrix A=r(table)
		foreach j of numlist ${b_age}/${e_age} {
			post combined_trans3 (2) (`j') (`i') (A[1,1]) (A[1,2]) (A[1,3]) (A[1,4]) (A[1,5])
		}

		* HR的自然对数beta及其标准误
		stcox i.HLS_5groups age_at_`i'_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e) nohr
		matrix u=r(table)
		post se ("female") ("HLS") (`i') (3) (u[1,1]) (u[1,2]) (u[1,3]) (u[1,4]) (u[1,5]) (u[2,1]) (u[2,2]) (u[2,3]) (u[2,4]) (u[2,5])

		* Table 1.
		loc cat_1=string(A[1,1],"%9.2f")+" ("+"Referent"+")"
		loc cat_2=string(A[1,2],"%9.2f")+" ("+string(A[5,2],"%9.2f")+"-"+string(A[6,2],"%9.2f")+")"
		loc cat_3=string(A[1,3],"%9.2f")+" ("+string(A[5,3],"%9.2f")+"-"+string(A[6,3],"%9.2f")+")"
		loc cat_4=string(A[1,4],"%9.2f")+" ("+string(A[5,4],"%9.2f")+"-"+string(A[6,4],"%9.2f")+")"
		loc cat_5=string(A[1,5],"%9.2f")+" ("+string(A[5,5],"%9.2f")+"-"+string(A[6,5],"%9.2f")+")"
		
		display "* Trend test"
		stcox c.HLS_5groups age_at_`i'_diagnosis i.highest_education i.marital_status i.heart_attack_FH i.stroke_FH i.cancer_FH i.menopause_inc_missing if is_female == 1, strata(region_code age_cox_strata) cformat(%9.2f) pformat(%9.2e)
		matrix B=r(table)
		loc trend=string(B[1,1],"%9.2f")+" ("+string(B[5,1],"%9.2f")+"-"+string(B[6,1],"%9.2f")+")"
		loc p=string(B[4,1],"%9.2e")
		post table1 ("Female") ("`i'") ("Disease → death") (`case') (`rate_1dec') ("`cat_1'") ("`cat_2'") ("`cat_3'") ("`cat_4'") ("`cat_5'") ("`trend'") ("`p'")
	}
	postclose combined_trans3
	postclose table1
	postclose se
	log close

	cd D:\SunQF\03-healthy_life_expectancy\Check\Cox
	use combined_trans3.dta, clear
	merge 1:1 sex age incidence using combined_trans1.dta, keep(match) nogenerate
	merge 1:1 sex age incidence using combined_trans2.dta, keep(match) nogenerate
	gen life = "HLS"
	save cox_combined_`datas'.dta, replace

	erase combined_trans1.dta
	erase combined_trans2.dta
	erase combined_trans3.dta

	* Table 1.
	use "D:\SunQF\03-healthy_life_expectancy\3 分析结果\Table-1_`datas'.dta", clear
	export excel using "D:\SunQF\03-healthy_life_expectancy\3 分析结果\Table-1_`datas'.xlsx", firstrow(variables) replace
	erase "D:\SunQF\03-healthy_life_expectancy\3 分析结果\Table-1_`datas'.dta"

	* For CI
	use "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_cox_combined_`datas'.dta", clear
	export excel using "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_cox_combined_`datas'.xlsx", firstrow(variables) replace
	erase "D:\SunQF\03-healthy_life_expectancy\Check\CI\se_cox_combined_`datas'.dta"
}


