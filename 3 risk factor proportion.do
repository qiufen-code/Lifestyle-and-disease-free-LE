*** Aim: sex- and age-specific risk factor prevalence for each transition

*** specify the database for use
loc prefix main

*** main code
foreach datas of loc prefix {

	* 不同分析中关注的疾病定义和生活方式因素
	if "`datas'" == "main" {
		loc life HLS_5groups smoking_5groups alcohol_5groups PA_5groups diet_5groups obesity_5groups
	}
	if "`datas'" != "main" {
		loc life HLS_5groups
	}


	postfile prop sex age str20 life incidence trans p1 p2 p3 p4 p5 using "D:\SunQF\03-healthy_life_expectancy\Check\proportion\follow_up\proportion_`datas'.dta", replace
	postfile prop_se sex age str20 life incidence trans p1 p2 p3 p4 p5 se1 se2 se3 se4 se5 using "D:\SunQF\03-healthy_life_expectancy\Check\CI\proportion_se_`datas'.dta", replace

	foreach group of loc life {

		* 生活方式组合关注多种疾病状态定义，单一生活方式只关注疾病组合定义
		if "`datas'" == "main" {
			loc incidence 99999 16418 16419 16420
		}
		else if "`datas'" == "T2D" {
			loc incidence 00000
		}
		else {
			loc incidence 99999
		}


		loc sublife = substr("`group'", 1, strlen("`group'")-8)

		foreach i of loc incidence {
			//* Transition 1
			cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
			use `datas'_`i'.dta, clear
			keep if _trans == 1
			stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
			stsplit age, every(5)
				replace age=80 if age>80                          //80岁之后作为一个age group
				
			* male
			sum age if is_female==0
			loc min = r(min)
			loc max = r(max)
			levelsof `group' if is_female==0, loc(level) clean    //单一生活方式及其组合的level都是1-5

			proportion `group' if is_female==0, over(age) cformat(%9.3f)
			mat A=r(table)
			foreach age of numlist `min'(5)`max' {
				foreach l of loc level {
					loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
					loc p`age'_`l'=A[1,`col']
					loc se`age'_`l'=A[2,`col']
				}
			}
			
			foreach age of numlist `min'(5)`max' {
				post prop (1) (`age') ("`sublife'") (`i') (1) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5')
				post prop_se (1) (`age') ("`sublife'") (`i') (1) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5') (`se`age'_1') (`se`age'_2') (`se`age'_3') (`se`age'_4') (`se`age'_5')
			}

			* female
			sum age if is_female==1
			loc min = r(min)
			loc max = r(max)
			levelsof `group' if is_female==1, loc(level) clean

			proportion `group' if is_female==1, over(age) cformat(%9.3f)
			mat A=r(table)
			foreach age of numlist `min'(5)`max' {
				foreach l of loc level {
					loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
					loc p`age'_`l'=A[1,`col']
					loc se`age'_`l'=A[2,`col']
				}
			}
			
			foreach age of numlist `min'(5)`max' {
				post prop (2) (`age') ("`sublife'") (`i') (1) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5')
				post prop_se (2) (`age') ("`sublife'") (`i') (1) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5') (`se`age'_1') (`se`age'_2') (`se`age'_3') (`se`age'_4') (`se`age'_5')
			}
			
			//* Transition 2
			cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
			use `datas'_`i'.dta, clear
			keep if _trans == 2
			stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
			stsplit age, every(5)
				replace age=80 if age>80

			* male
			sum age if is_female==0
			loc min = r(min)
			loc max = r(max)
			levelsof `group' if is_female==0, loc(level) clean

			proportion `group' if is_female==0, over(age) cformat(%9.3f)
			mat A=r(table)
			foreach age of numlist `min'(5)`max' {
				foreach l of loc level {
					loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
					loc p`age'_`l'=A[1,`col']
					loc se`age'_`l'=A[2,`col']
				}
			}
			
			foreach age of numlist `min'(5)`max' {
				post prop (1) (`age') ("`sublife'") (`i') (2) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5')
				post prop_se (1) (`age') ("`sublife'") (`i') (2) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5') (`se`age'_1') (`se`age'_2') (`se`age'_3') (`se`age'_4') (`se`age'_5')
			}
			
			* female
			sum age if is_female==1
			loc min = r(min)
			loc max = r(max)
			levelsof `group' if is_female==1, loc(level) clean

			proportion `group' if is_female==1, over(age) cformat(%9.3f)
			mat A=r(table)
			foreach age of numlist `min'(5)`max' {
				foreach l of loc level {
					loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
					loc p`age'_`l'=A[1,`col']
					loc se`age'_`l'=A[2,`col']
				}
			}
			
			foreach age of numlist `min'(5)`max' {
				post prop (2) (`age') ("`sublife'") (`i') (2) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5')
				post prop_se (2) (`age') ("`sublife'") (`i') (2) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5') (`se`age'_1') (`se`age'_2') (`se`age'_3') (`se`age'_4') (`se`age'_5')
			}
			
			* Transition 3
			cd D:\SunQF\03-healthy_life_expectancy\Check\mstate
			use `datas'_`i'.dta, clear
			keep if _trans == 3
			stset _stop, id(csid) enter(_start) origin(time dob_anon) failure(_status==1) scale(365.25)
			stsplit age, every(5)
				replace age=floor((_start-dob_anon)/365.25/5)*5 if age==.
				replace age=80 if age>80

			* male
			sum age if is_female==0
			loc min = r(min)
			loc max = r(max)
			levelsof `group' if is_female==0, loc(level) clean

			proportion `group' if is_female==0, over(age) cformat(%9.3f)
			mat A=r(table)
			foreach age of numlist `min'(5)`max' {
				foreach l of loc level {
					loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
					loc p`age'_`l'=A[1,`col']
					loc se`age'_`l'=A[2,`col']
				}
			}
			
			foreach age of numlist `min'(5)`max' {
				post prop (1) (`age') ("`sublife'") (`i') (3) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5')
				post prop_se (1) (`age') ("`sublife'") (`i') (3) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5') (`se`age'_1') (`se`age'_2') (`se`age'_3') (`se`age'_4') (`se`age'_5')
			}
			
			* female
			sum age if is_female==1
			loc min = r(min)
			loc max = r(max)
			levelsof `group' if is_female==1, loc(level) clean

			proportion `group' if is_female==1, over(age) cformat(%9.3f)
			mat A=r(table)
			foreach age of numlist `min'(5)`max' {
				foreach l of loc level {
					loc col=(`age'-`min')/5 + 1 + (`l'-1)*((`max'-`min')/5+1)
					loc p`age'_`l'=A[1,`col']
					loc se`age'_`l'=A[2,`col']
				}
			}
			
			foreach age of numlist `min'(5)`max' {
				post prop (2) (`age') ("`sublife'") (`i') (3) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5')
				post prop_se (2) (`age') ("`sublife'") (`i') (3) (`p`age'_1') (`p`age'_2') (`p`age'_3') (`p`age'_4') (`p`age'_5') (`se`age'_1') (`se`age'_2') (`se`age'_3') (`se`age'_4') (`se`age'_5')
			}
		}
	}
	postclose prop
	postclose prop_se


	* 对proportion数据集进行转置
	use "D:\SunQF\03-healthy_life_expectancy\Check\proportion\follow_up\proportion_`datas'.dta", clear
	reshape wide p@1 p@2 p@3 p@4 p@5, i(sex age life incidence) j(trans)
	sort sex age life incidence
	save "D:\SunQF\03-healthy_life_expectancy\Check\proportion\follow_up\proportion_`datas'.dta", replace

	* For CI
	use "D:\SunQF\03-healthy_life_expectancy\Check\CI\proportion_se_`datas'.dta", clear
	export excel using "D:\SunQF\03-healthy_life_expectancy\Check\CI\proportion_se_`datas'.xlsx", firstrow(variables) replace
	erase "D:\SunQF\03-healthy_life_expectancy\Check\CI\proportion_se_`datas'.dta"
}

