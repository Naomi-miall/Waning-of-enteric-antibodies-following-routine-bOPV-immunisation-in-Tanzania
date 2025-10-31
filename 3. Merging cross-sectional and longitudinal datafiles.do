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
tab pv1_stool_positive visit if sampletype==1 & flag_ageoutlier==0, col
bysort visit: summ logPV1_stool if sampletype==1 & flag_ageoutlier==0, d

tab pv3_stool_positive visit if sampletype==1 & flag_ageoutlier==0, col
bysort visit: summ logPV3_stool if sampletype==1 & flag_ageoutlier==0, d


bysort visit: summ logPV1_serum if sampletype==1 & flag_ageoutlier==0, d
bysort visit: summ logPV3_serum if sampletype==1 & flag_ageoutlier==0, d


tab pv3_stool_responses visit if flag_ageoutlier==0 & sampletype==1, col
tab pv3_stool_responses visit if flag_ageoutlier==0 & sampletype==2, col
tab pv3_stool_responses visit if flag_ageoutlier==0, col

tab pv1_stool_responses visit if flag_ageoutlier==0 & sampletype==1, col
tab pv1_stool_responses visit if flag_ageoutlier==0 & sampletype==2, col
tab pv1_stool_responses visit if flag_ageoutlier==0, col

///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
*Regressions: Estimating half lives:

*Linear regression for PV1 between 1 month and 6 months
use "merged_dataset_long", clear
keep if visit>1 & flag_ageoutlier==0
mixed logPV1_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid: //coefficient is -.0096988
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days //coefficient is -.0096988
summ logPV1_stool if visit ==2, d //median is 2.424242 so half value is 1.212121
display 1.212121/.0096988 //half life = 125.0 days

//Including the cross-sectional sample also
use "Both_datasets_long", clear
keep if visit>1 & visit<5 & flag_ageoutlier==0
mixed logPV1_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid:
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days //coefficient is  -.0085907
summ logPV1_stool if visit ==2, d //median is 1.535354 so half value is 0.767677
display 0.767677/.0085907 //half life = 89.4 days

use "merged_dataset_long", clear
keep if visit>1 & flag_ageoutlier==0
mixed logPV3_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV3_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV3_stool if visit ==2, d
display ( 2.820707)/.0196951 

//Including the cross-sectional sample also
use "Both_datasets_long", clear
keep if visit>1 & visit<5 & flag_ageoutlier==0
mixed logPV3_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV3_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV3_stool if visit ==2, d 
display 2.608586/.0179969 //half life = 145


//High baseline
use "merged_dataset_long", clear
keep if visit>1 & flag_ageoutlier==0 & high_neut_baseline_pv1==1
mixed logPV1_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV1_stool if visit ==2, d
display (3.025253/2)/.0110366   

use "merged_dataset_long", clear
keep if visit>1 & flag_ageoutlier==0 & high_neut_baseline_pv3==1
mixed logPV3_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV3_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV3_stool if visit ==2, d
display ( 5.888889/2)/.0216188


//Low baseline 
use "merged_dataset_long", clear
keep if visit>1 & flag_ageoutlier==0 & high_neut_baseline_pv1==0
mixed logPV1_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV1_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV1_stool if visit ==2, d
display (1.606061/2)/.0091918


use "merged_dataset_long", clear
keep if visit>1 & flag_ageoutlier==0 & high_neut_baseline_pv3==0
mixed logPV3_stool follow_up_months femalesex age_opv14_weeks || pid: follow_up_months
mixed logPV3_stool follow_up_days femalesex age_opv14_weeks || pid: follow_up_days
summ logPV3_stool if visit ==2, d
display (3.414141/2)/.0146626 

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

keep if (sampletype==2 & visit==1 & flag_ageoutlier==0) | (sampletype==1 & visit<5 & flag_ageoutlier==0)

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
keep if sampletype==2 & visit==1 & flag_ageoutlier==0

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

keep if (sampletype==2 & visit==1 & flag_ageoutlier==0) | (sampletype==1 & visit<5 & flag_ageoutlier==0)
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

keep if (sampletype==2 & visit==1 & flag_ageoutlier==0) | (sampletype==1 & visit<5 & flag_ageoutlier==0)
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

////////////////////////////////////////////////////////////
*correlations
use "Both_datasets_long", clear
*using month 1 for longitudinal:
keep if sampletype==1 | visit==2
tab visit sampletype
save "correlations_dataset", replace

*STOOL AND SERUM:
//NEED TO THINK ABOUT MANAGING REPEAT OBS FOR LONG SAMPLE
spearman logPV1_stool logPV1_serum, exact //rho = 0.33, p=0.001
spearman logPV3_stool logPV3_serum, exact //rho = 0.51, p<0.001

spearman logPV1_stool logPV3_stool 
spearman logPV1_serum logPV3_serum

///////////////////////////////////////////////////////////////
*Creating a dataset for the any_infections 
use "Both_datasets_long", clear
keep if visit==2 & any_infection!=. & flag_ageoutlier==0

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
