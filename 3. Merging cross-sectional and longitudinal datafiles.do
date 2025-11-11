/* Duration of polio antibody responses following bOPV vaccination in Tanzania.
  Alexandra Codeiro and Naomi Miall
  September 2025
*/
clear
set more off

*Filepath for Naomi:
cd "C:\Users\NaomiMiall\Documents\Polio work\Tanzania"
///////////////////////////////////////////////////////////////

*Append cross-sectional
use "merged_dataset_long", clear
append using cross-sectional_clean

foreach x in 1 2 3 {
	replace  pv`x'_serum = pv`x'neuttiterserum if sampletype==2
	drop pv`x'neuttiterserum
}

gen follow_up_weeks_integer = int(follow_up_days/7)
gen follow_up_months_integer = int(follow_up_days/30)

//Those followed-up > 7 days after bOPV4 are equally likley to be PV1/PV3 positive as those followed up < 7 days after.
tab pv3_stool_positive follow_up_weeks_integer if visit==1, col
tab pv1_stool_positive follow_up_weeks_integer if visit==1, col

tab follow_up_months_integer visit, m
tab follow_up_weeks_integer visit, m
//Weeks for long are => 0-1 for baseline; 2-6 for month 1; 11-18 for month 3; 24-31 for month 6

gen age_bopv14w = int((bopv14w_date - infant_dob)/7)
tab age_bopv14w
hist age_bopv14w, bin(100)

*gen visit_cross = .
*replace visit_cross = followup_weeks if sampletype==1
*recode visit_cross 0/6=1 7/21=2 22/32=3 46/65=4 96/111=5 140/170=6 400/550=7

gen visit_cross = .
replace visit_cross = follow_up_months_integer if sampletype==1
recode visit_cross 0/1=2 2/4=3 5/7=4 10/15=5 22/26=6 34/38=7 97/119=8
label define visit_cross 2 "Month 1" 3 "Month 3" 4 "Month 6" 5 "Year 1" 6 "Year 2" 7 "Year 3" 8 "Year 8-9", replace
label values visit_cross visit_cross
tab visit_cross
tab follow_up_weeks_integer visit_cross if sampletype==1

replace visit = visit_cross if sampletype==1

*variable for average 'day' for each visit (needed for summary stats on plots)
drop summ_visit_day
bysort visit: summ follow_up_days
gen summ_visit_day = .
replace summ_visit_day = 0 if visit==1 
replace summ_visit_day = 30 if visit==2
replace summ_visit_day = 90 if visit==3
replace summ_visit_day = 180 if visit==4
replace summ_visit_day = 365 if visit==5
replace summ_visit_day = 730 if visit==6
replace summ_visit_day = 1095 if visit==7
replace summ_visit_day = 3103 if visit==8

replace pv1_stool_responses = pv1_stool_positive if sampletype==1
replace pv1_stool_responses_exclbl = pv1_stool_positive if sampletype==1
replace pv3_stool_responses = pv3_stool_positive if sampletype==1
replace pv3_stool_responses_exclbl = pv3_stool_positive if sampletype==1

save "Both_datasets_long", replace
use "Both_datasets_long", clear

///////////////////////////////////////////////////////////////
tab pv1_stool_positive visit if sampletype==1, col
bysort visit: summ logPV1_stool if sampletype==1, d

tab pv3_stool_positive visit if sampletype==1, col
bysort visit: summ logPV3_stool if sampletype==1, d


bysort visit: summ logPV1_serum if sampletype==1, d
bysort visit: summ logPV3_serum if sampletype==1, d


tab pv3_stool_responses visit if sampletype==1, col
tab pv3_stool_responses visit if sampletype==2, col
tab pv3_stool_responses visit, col

tab pv1_stool_responses visit if sampletype==1, col
tab pv1_stool_responses visit if sampletype==2, col
tab pv1_stool_responses visit if flag_ageoutlier==0, col

bysort visit: summ logPV1_IgA_stool if sampletype==1, d
bysort visit: summ logPV3_IgA_stool if sampletype==1, d
bysort visit: summ logPV1_IgA_serum if sampletype==1, d
bysort visit: summ logPV3_IgA_serum if sampletype==1, d
bysort visit: summ logPV1_IgG_serum if sampletype==1, d
bysort visit: summ logPV3_IgG_serum if sampletype==1, d

*Supplementary table 3
use "Both_datasets_long", clear
 keep if sampletype == 1

postfile handle str30 variable str20 visit_group N median ci_low ci_high using results_temp.dta, replace

local vars logPV1_IgA_stool logPV1_IgA_stool logPV1_IgA_stool logPV1_IgA_stool logPV1_IgA_stool logPV1_IgA_stool logPV1_IgA_stool logPV3_IgA_stool logPV3_IgA_stool logPV3_IgA_stool logPV3_IgA_stool logPV3_IgA_stool logPV3_IgA_stool logPV3_IgA_stool logPV1_IgA_serum logPV1_IgA_serum logPV1_IgA_serum logPV1_IgA_serum logPV1_IgA_serum logPV1_IgA_serum logPV1_IgA_serum logPV3_IgA_serum logPV3_IgA_serum logPV3_IgA_serum logPV3_IgA_serum logPV3_IgA_serum logPV3_IgA_serum logPV3_IgA_serum logPV1_IgG_serum logPV1_IgG_serum logPV1_IgG_serum logPV1_IgG_serum logPV1_IgG_serum logPV1_IgG_serum logPV1_IgG_serum logPV3_IgG_serum logPV3_IgG_serum logPV3_IgG_serum logPV3_IgG_serum logPV3_IgG_serum logPV3_IgG_serum logPV3_IgG_serum 
local visit 2 3 4 5 6 7 8 2 3 4 5 6 7 8 2 3 4 5 6 7 8 2 3 4 5 6 7 8 2 3 4 5 6 7 8 2 3 4 5 6 7 8
forvalues x=1/42 {
	local vars2: word `x' of `vars'
	local visit2: word `x' of `visit'
summarize `vars2' if visit==`visit2', detail
post handle ("`vars2'") ("`visit2'") (r(N)) (r(p50)) (r(p25)) (r(p75))
}
postclose handle
use results_temp.dta, clear
export excel using "supplementarytab3.xlsx", sheet("Sheet 1") replace

summarize logPV1_IgA_stool if visit==2, d
local N = r(N)
local p25 = r(p25)
display `p25'


use "Both_datasets_long", clear
///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
*Regressions: Estimating half lives:

*Linear regression for PV1 between 1 month and 6 months
use "both_datasets_long", clear
keep if sampletype==1 & visit<5 
regress logPV1_stool follow_up_months femalesex age_opv14_weeks //-0.1638856 
summ logPV1_stool if visit ==2, d //median is 1 - lowest possible already

use "both_datasets_long", clear
keep if sampletype==1 & visit<5 
regress logPV3_stool follow_up_months femalesex age_opv14_weeks //-.4120233
regress logPV3_stool follow_up_days femalesex age_opv14_weeks //.0137341
summ logPV3_stool if visit ==2, d //median is  3.555556  
display  (3.555556/2)/.0137341


//High baseline
use "merged_dataset_long", clear
keep if visit>1 & high_neut_baseline_pv1==1
mixed logPV1_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV1_stool if visit ==2, d
display (3.025253/2)/.0110366   

use "merged_dataset_long", clear
keep if visit>1 & high_neut_baseline_pv3==1
mixed logPV3_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV3_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV3_stool if visit ==2, d
display ( 5.959596 /2)/.0218261


//Low baseline 
use "merged_dataset_long", clear
keep if visit>1 & high_neut_baseline_pv1==0
mixed logPV1_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV1_stool if visit ==2, d
display (2.212121/2)/.0093252 


use "merged_dataset_long", clear
keep if visit>1 & high_neut_baseline_pv3==0
mixed logPV3_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV3_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV3_stool if visit ==2, d
display (2.777778/2)/.013512

////////////////////////////////////////////////////////////////////
/*Regressions with covariates:
Adjusting for: 
Set 1: age, sex
Set 2: age, sex, breastfeeding status, primary water source, and stunting.
Set 3: age, sex, anyinfection (includes schisto, parasites, biofire but not HIV exposure)
Set 4: age, sex, and each coinfection separately */ 
*for now adjusted for age as a binary but in future I'll adjust for age as continuous.

use "Both_datasets_long", clear
sort pid visit
list pid sampletype visit pv1_stool logPV1_stool pv1_stool_positive pv1_stool_responses pv1_stool_responses_exclbl in 1/50, sepby(pid)

keep if (sampletype==2 & visit==1) | (sampletype==1 & visit<5)

logistic pv1_stool_responses_exclbl bf_4months ageweek femalesex, or
logistic pv1_stool_responses_exclbl watersource_bin ageweek femalesex, or
logistic pv1_stool_responses_exclbl stunting ageweek femalesex, or
logistic pv1_stool_responses_exclbl any_infection ageweek femalesex, or
logistic pv1_stool_responses_exclbl hiv_exposure ageweek femalesex, or
logistic pv1_stool_responses_exclbl coinfect_caa ageweek femalesex, or

logistic pv1_stool_responses_exclbl coinfect_giardia ageweek femalesex, or
logistic pv1_stool_responses_exclbl coinfect_enteroinvasive ageweek femalesex, or
logistic pv1_stool_responses_exclbl coinfect_otherbacteria ageweek femalesex, or
logistic pv1_stool_responses_exclbl ageweek femalesex coinfect_virus, or


logistic pv3_stool_responses_exclbl bf_4months ageweek femalesex, or
logistic pv3_stool_responses_exclbl watersource_bin ageweek femalesex, or
logistic pv3_stool_responses_exclbl stunting ageweek femalesex, or
logistic pv3_stool_responses_exclbl any_infection ageweek femalesex, or
logistic pv3_stool_responses_exclbl hiv_exposure ageweek femalesex, or
logistic pv3_stool_responses_exclbl coinfect_caa ageweek femalesex, or

logistic pv3_stool_responses_exclbl coinfect_giardia ageweek femalesex, or
logistic pv3_stool_responses_exclbl coinfect_enteroinvasive ageweek femalesex, or
logistic pv3_stool_responses_exclbl coinfect_otherbacteria ageweek femalesex, or
logistic pv3_stool_responses_exclbl ageweek femalesex coinfect_virus, or

use "Both_datasets_long", clear
drop if visit!=2

foreach x in bf_4months watersource_bin stunting any_infection coinfect* {
	logistic pv1_stool_positive `x' ageweek femalesex, or
	logistic pv3_stool_positive `x' ageweek femalesex, or
}

//Longitudinal only - adjusting for baseline
use "Both_datasets_long", clear
keep if sampletype==2 & visit==1

postfile handle str20 outcome coef pval ci_low ci_high using results_temp.dta, replace

foreach x in bf_4months watersource_bin stunting hiv_exposure any_infection coinfect_caa coinfect_giardia coinfect_enteroinvasive coinfect_otherbacteria coinfect_virus {
	logistic pv1_stool_responses_exclbl `x' ageweek femalesex high_neut_baseline_pv1, or
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle ("pv1`x'") (`b') (`p') (`ci_l') (`ci_h')
	
	logistic pv3_stool_responses_exclbl `x' ageweek femalesex high_neut_baseline_pv3, or
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle ("pv3`x'") (`b') (`p') (`ci_l') (`ci_h')
}
postclose handle
use results_temp.dta, clear
export excel using "reg_long_results2.xlsx", replace

//Looking at detectable 
use "Both_datasets_long", clear
gen detectable_pv1_prep = pv1_stool_detectable
replace detectable_pv1_prep=. if visit==1
bysort pid: egen detectable_pv1 = max(detectable_pv1_prep)
list pid visit sampletype detectable_pv1 pv1_stool_detectable in 1/50, sepby(pid)

gen detectable_pv3_prep = pv3_stool_detectable
replace detectable_pv3_prep=. if visit==1
bysort pid: egen detectable_pv3 = max(detectable_pv3_prep)
list pid visit sampletype detectable_pv3 pv3_stool_detectable in 1/50, sepby(pid)

keep if (sampletype==2 & visit==1 ) | (sampletype==1 & visit<5)
postfile handle2 str20 outcome str30 predictor coef pval ci_low ci_high using results_temp2.dta, replace
foreach x in bf_4months watersource_bin stunting hiv_exposure any_infection coinfect_caa coinfect_enteroinvasive coinfect_otherbacteria coinfect_virus {
	logistic detectable_pv1  `x' ageweek femalesex, or
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle2 ("pv1") ("`x'") (`b') (`p') (`ci_l') (`ci_h')
	
	logistic detectable_pv3  `x' ageweek femalesex, or
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle2 ("pv3") ("`x'") (`b') (`p') (`ci_l') (`ci_h')
}
postclose handle2
use results_temp2.dta, clear
export excel using "reg_long_results3.xlsx", replace


//Looking at max neut titre measured after baseline
use "Both_datasets_long", clear
foreach x in 1 2 3 {
	gen logpv`x'_stool_exclbl = logPV`x'_stool
	replace logpv`x'_stool_exclbl = 0 if visit==1
	bysort pid: egen max_logpv`x'_stool_exclbl = max(logpv`x'_stool_exclbl)
	drop logpv`x'_stool_exclbl
}

sort pid visit
list pid visit sampletype logPV1_stool max_logpv1_stool_exclbl logPV3_stool max_logpv3_stool_exclbl in 1/200 if sampletype==2, sepby(pid)
list pid visit sampletype logPV1_stool max_logpv1_stool_exclbl logPV3_stool max_logpv3_stool_exclbl in 1/50, sepby(pid)

keep if (sampletype==2 & visit==1) | (sampletype==1 & visit<5)
postfile handle3 str20 outcome str30 predictor coef pval ci_low ci_high using results_temp3.dta, replace
foreach x in bf_4months watersource_bin stunting hiv_exposure any_infection coinfect_caa coinfect_enteroinvasive coinfect_otherbacteria coinfect_virus {
	regress max_logpv1_stool_exclbl  `x' ageweek femalesex
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle3 ("pv1") ("`x'") (`b') (`ci_l') (`ci_h') (`p')
	
	regress max_logpv3_stool_exclbl  `x' ageweek femalesex
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle3 ("pv3") ("`x'") (`b') (`ci_l') (`ci_h') (`p')
}
postclose handle3
use results_temp3.dta, clear
export excel using "reg_long_results4.xlsx", replace



//Looking at neut titre at 1 month
use "Both_datasets_long", clear
sort pid visit
keep if visit==2
postfile handle4 str20 outcome str30 predictor coef pval ci_low ci_high using results_temp4.dta, replace
foreach x in bf_4months watersource_bin stunting hiv_exposure any_infection coinfect_caa coinfect_enteroinvasive coinfect_otherbacteria coinfect_virus {
	logistic pv1_stool_positive  `x' ageweek femalesex, or
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle4 ("pv1") ("`x'") (`b') (`ci_l') (`ci_h') (`p')
	
	logistic pv3_stool_positive  `x' ageweek femalesex, or
	local b = _b[`x']
    local ci_l = _b[`x'] - 1.96 * _se[`x']
    local ci_h = _b[`x'] + 1.96 * _se[`x']
	testparm `x'
	local p = r(p)
	post handle4 ("pv3") ("`x'") (`b') (`ci_l') (`ci_h') (`p')
}
postclose handle4
use results_temp4.dta, clear
export excel using "reg_long_results5.xlsx", replace

////////////////////////////////////////////////////////////
*correlations
use "Both_datasets_long", clear
*using month 1 for longitudinal:
keep if sampletype==1 | visit==2
tab visit sampletype
save "correlations_dataset", replace

/**STOOL AND SERUM:
//NEED TO THINK ABOUT MANAGING REPEAT OBS FOR LONG SAMPLE
spearman logPV1_stool logPV1_serum, exact //rho = 0.33, p=0.001
spearman logPV3_stool logPV3_serum, exact //rho = 0.51, p<0.001

spearman logPV1_stool logPV3_stool 
spearman logPV1_serum logPV3_serum
*/


///////////////////////////////////////////////////////////////
*Creating a dataset for the any_infections 
use "Both_datasets_long", clear
keep if visit==2 & any_infection!=.

rename logPV1_stool logPV_1
rename logPV3_stool logPV_2
rename logPV1_serum logPV_3
rename logPV3_serum logPV_4

reshape long logPV_, i(pid) j(assay)
sort pid assay
label define assay 1 "PV1 stool" 2"PV3 stool" 3"PV1 serum" 4"PV3 serum"
label values assay assay
list pid visit sampletype assay logPV_ in 1/20, sepby(pid)

save "infections_dataset", replace



///////////////////////////////////////////////////////////////
*Creating a dataset for the box_plots 
use "Both_datasets_long", clear
keep pid visit sampletype logPV*

rename logPV1_stool logPV1
rename logPV3_stool logPV2
rename logPV1_serum logPV3
rename logPV3_serum logPV4
rename logPV1_IgA_stool logPV5
rename logPV3_IgA_stool logPV6
rename logPV1_IgA_serum logPV7
rename logPV3_IgA_serum logPV8
rename logPV1_IgG_serum logPV9
rename logPV3_IgG_serum logPV10

reshape long logPV, i(pid visit) j(assay)
sort pid assay
label define assay 1 "PV1 stool neut" 2"PV3 stool neut" 3"PV1 serum neut" 4"PV3 serum neut" 5"PV1 stool IgA" 6"PV3 stool IgA" 7"PV1 serum IgA" 8"PV3 serum IgA" 9"PV1 serum IgG" 10"PV3 serum IgG"
label values assay assay
list pid visit sampletype assay logPV in 1/20, sepby(pid)

save "box_plots", replace



//Box plot file for the PV2 box plots
use "Both_datasets_long", clear
keep pid visit sampletype logPV*

rename logPV2_stool logPV1
rename logPV2_serum logPV2
rename logPV2_IgA_stool logPV3
rename logPV2_IgA_serum logPV4
rename logPV2_IgG_serum logPV5

reshape long logPV, i(pid visit) j(assay)
sort pid assay
label define assay 1 "PV2 stool neut" 2 "PV2 serum neut" 3 "PV2 stool IgA" 4 "PV2 serum IgA" 5 "PV2 serum IgG"
label values assay assay
list pid visit sampletype assay logPV in 1/20, sepby(pid)

save "box_plots_PV2", replace
