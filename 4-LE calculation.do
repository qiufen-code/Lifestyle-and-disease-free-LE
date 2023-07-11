* Disease-free life expectancy

foreach s of numlist 1/5 {
    * male
    use D:\SunQF\03-healthy_life_expectancy\Check\all\all99999_HLS_main.dta, clear
    keep if sex == 1
    drop if age < 40 | age > 90
    keep age hazard1`s' hazard2`s' hazard3`s'
    rename hazard1`s' hazard12`s'
    rename hazard2`s' hazard13`s'
    rename hazard3`s' hazard23`s'
    gen hazard21`s' = 0
    order age hazard12`s' hazard13`s' hazard21`s' hazard23`s'
    mslt, l0(100000 0 0) death censored

    * female
    use D:\SunQF\03-healthy_life_expectancy\Check\all\all99999_HLS_main.dta, clear
    keep if sex == 2
    drop if age < 40 | age > 90
    keep age hazard1`s' hazard2`s' hazard3`s'
    rename hazard1`s' hazard12`s'
    rename hazard2`s' hazard13`s'
    rename hazard3`s' hazard23`s'
    gen hazard21`s' = 0
    order age hazard12`s' hazard13`s' hazard21`s' hazard23`s'
    mslt, l0(100000 0 0) death censored
}
