* Disease-free life expectancy
* Notation for installation of "moremata" in stata：
* 1、下载moremata.zip
* 2、解压在c:\temp (注意是将moremata中的所有子文件解压到temp文件夹)
* 3、net from "c:\temp"
* 4、net install moremata, replace

*** specify the database for use
loc prefix main

*** main code
foreach datas of loc prefix {

  * 不同分析中关注的疾病定义和生活方式因素
  if "`datas'" == "main" {
    loc life HLS_5groups smoking_5groups alcohol_5groups PA_5groups diet_5groups obesity_5groups
    loc mid
  }
  else if "`datas'" != "main" {
    loc life HLS_5groups
    loc mid combined_
  }

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

    foreach i of loc incidence {
      loc sublife = substr("`group'", 1, strlen("`group'")-8)

      * 将rate和proportion数据集合并
      use "D:\SunQF\03-healthy_life_expectancy\Check\proportion\follow_up\proportion_`datas'", clear
      keep if incidence == `i' & life == "`sublife'"
      drop incidence life
      tempfile proportion
      save "`proportion'"

      use D:\SunQF\03-healthy_life_expectancy\Check\trans_rate\rate`i'_`datas'.dta, clear
      merge 1:1 sex age using "`proportion'", nogenerate
      gen incidence = `i'
      gen life = "`sublife'"

      foreach p of numlist 1/3 { 
        foreach l of numlist 1/5 {
          replace p`p'`l' = . if age == 90
          replace p`p'`l' = p`p'`l'[_n-1] if _n > 1 & p`p'`l' == .     //填充90岁以上的proportion
        }
      }
      preserve

      * 与cox数据集合并
      use D:\SunQF\03-healthy_life_expectancy\Check\Cox\cox_`mid'`datas'.dta, clear
      keep if incidence == `i' & life == "`sublife'"
      tempfile cox
      save "`cox'", replace

      * 估算每一种生活方式水平(level)者的transition rates
      restore
      merge 1:1 sex age incidence using "`cox'", keep(match) nogenerate
      foreach trans of numlist 1/3 {
        foreach l of numlist 1/5 {                                       //在根据生活方式因素其一做亚组分析时，需要注意修改水平数
          tempvar sum`trans'`l'
          gen `sum`trans'`l'' = p`trans'`l'*hr`trans'`l'
        }
        tempvar sum`trans'
        egen `sum`trans'' = rowtotal(`sum`trans'1'-`sum`trans'5')
      }

      foreach trans of numlist 1/3 {
        foreach l of numlist 1/5 {
          gen hazard`trans'`l' = hr`trans'`l' * (hazard`trans'/`sum`trans'')
        }
      }

      cd "D:\SunQF\03-healthy_life_expectancy\Check\all"
      save all`i'_`sublife'_`datas'.dta, replace
    } 
  }
}


*******************************************************************************************
************************************                   ************************************
************************************ Health expectancy ************************************
************************************                   ************************************
*******************************************************************************************

*** main code for estimation of HLE
loc prefix main
global b_age = 40
global e_age = 90
set graphics off

foreach datas of loc prefix {

  * 敏感性分析和亚组分析中只关注生活方式组合
  if "`datas'" == "main" & ${b_age} == 40 {
    loc life HLS smoking alcohol PA diet obesity
  }
  else if "`datas'" != "main" {
    loc life HLS
  }

  cap postclose LE
  postfile LE str20 (sex exposure group disease status years) using "D:\SunQF\03-healthy_life_expectancy\Check\figure\LE_at${b_age}_`datas'.dta", replace

  foreach sublife of loc life {

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

    foreach i of loc incidence {

      foreach s of numlist 1/5 {
        * male
        use D:\SunQF\03-healthy_life_expectancy\Check\all\all`i'_`sublife'_`datas'.dta, clear
        keep if sex == 1
        drop if age < ${b_age} | age > ${e_age}
        keep age hazard1`s' hazard2`s' hazard3`s'
        rename hazard1`s' hazard12`s'
        rename hazard2`s' hazard13`s'
        rename hazard3`s' hazard23`s'
        gen hazard21`s' = 0
        order age hazard12`s' hazard13`s' hazard21`s' hazard23`s'
        mslt, l0(100000 0 0) death censored                                         //采用默认的linear法估计概率，结尾年龄区间采用闭区间
        
        mat A=ei_x
        loc healthy`s'=string(A[1,1],"%9.1f")                                       //${b_age}时的不同状态(healthy、disease)下的期望寿命
        loc disease`s'=string(A[1,2],"%9.1f")
        post LE ("Male") ("`sublife'") ("`s'") ("`i'") ("healthy") ("`healthy`s''")
        post LE ("Male") ("`sublife'") ("`s'") ("`i'") ("disease") ("`disease`s''")
        
        if "`datas'" == "main" & "`sublife'" == "HLS" & `i' == 99999 {              //${b_age}至${e_age}期间时的全部期望寿命水平
          foreach age of numlist ${b_age}/${e_age} {
            loc row = `age' - ${b_age} + 1
            loc Mhealthy`age'_`s'=string(A[`row',1],"%9.1f")
            loc Mdisease`age'_`s'=string(A[`row',2],"%9.1f")
          }
        }
        
        * female
        use D:\SunQF\03-healthy_life_expectancy\Check\all\all`i'_`sublife'_`datas'.dta, clear
        keep if sex == 2
        drop if age < ${b_age} | age > ${e_age}
        keep age hazard1`s' hazard2`s' hazard3`s'
        rename hazard1`s' hazard12`s'
        rename hazard2`s' hazard13`s'
        rename hazard3`s' hazard23`s'
        gen hazard21`s' = 0
        order age hazard12`s' hazard13`s' hazard21`s' hazard23`s'
        mslt, l0(100000 0 0) death censored
        
        mat A=ei_x
        loc healthy`s'=string(A[1,1],"%9.1f")
        loc disease`s'=string(A[1,2],"%9.1f")
        post LE ("Female") ("`sublife'") ("`s'") ("`i'") ("healthy") ("`healthy`s''")
        post LE ("Female") ("`sublife'") ("`s'") ("`i'") ("disease") ("`disease`s''")
        
        if "`datas'" == "main" & "`sublife'" == "HLS" & `i' == 99999 {
          foreach age of numlist ${b_age}/${e_age} {
            loc row = `age' - ${b_age} + 1
            loc Fhealthy`age'_`s'=string(A[`row',1],"%9.1f")
            loc Fdisease`age'_`s'=string(A[`row',2],"%9.1f")
          }
        }
      }
        
      * figure s3 
      if "`datas'" == "main" & "`sublife'" == "HLS" & `i' == 99999 {
        cap postclose all_ages
        postfile all_ages str20 sex age str20 status ex1 ex2 ex3 ex4 ex5 using "D:\SunQF\03-healthy_life_expectancy\Check\figure\all_ages_`datas'.dta", replace
        
        foreach age of numlist ${b_age}/${e_age} {
          post all_ages ("Male") (`age') ("healthy") (`Mhealthy`age'_1') (`Mhealthy`age'_2') (`Mhealthy`age'_3') (`Mhealthy`age'_4') (`Mhealthy`age'_5')
          post all_ages ("Male") (`age') ("disease") (`Mdisease`age'_1') (`Mdisease`age'_2') (`Mdisease`age'_3') (`Mdisease`age'_4') (`Mdisease`age'_5')
        }
        
        foreach age of numlist ${b_age}/${e_age} {
          post all_ages ("Female") (`age') ("healthy") (`Fhealthy`age'_1') (`Fhealthy`age'_2') (`Fhealthy`age'_3') (`Fhealthy`age'_4') (`Fhealthy`age'_5')
          post all_ages ("Female") (`age') ("disease") (`Fdisease`age'_1') (`Fdisease`age'_2') (`Fdisease`age'_3') (`Fdisease`age'_4') (`Fdisease`age'_5')
        }
        
        postclose all_ages
        use "D:\SunQF\03-healthy_life_expectancy\Check\figure\all_ages_`datas'.dta", clear
        export excel using "D:\SunQF\03-healthy_life_expectancy\Check\figure\all_ages_`datas'.xlsx", firstrow(variables) replace
        erase "D:\SunQF\03-healthy_life_expectancy\Check\figure\all_ages_`datas'.dta"
      }
    }
  }
  postclose LE
        
  * 结果输出
  use "D:\SunQF\03-healthy_life_expectancy\Check\figure\LE_at${b_age}_`datas'.dta", clear
  export excel using "D:\SunQF\03-healthy_life_expectancy\Check\figure\LE_at${b_age}_`datas'.xlsx", firstrow(variables) replace
  erase "D:\SunQF\03-healthy_life_expectancy\Check\figure\LE_at${b_age}_`datas'.dta"
}




