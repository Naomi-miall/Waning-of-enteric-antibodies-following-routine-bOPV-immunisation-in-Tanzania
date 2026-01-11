/* Duration of polio antibody responses following bOPV vaccination in Tanzania.
  Alexandra Codeiro and Naomi Miall
  September 2025
*/
clear
set more off

*Filepath for Naomi:
cd "C:\Users\NaomiMiall\Documents\Polio work\Tanzania"
///////////////////////////////////////////////////////////////
use "merged_dataset_long", clear
append using "cross-sectional_clean"
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

*Generate visit number for cross-sectional participants
gen visit_cross = .
replace visit_cross = follow_up_months_integer if sampletype==1
recode visit_cross 0/1=2 2/4=3 5/7=4 10/15=5 22/26=6 34/38=7 97/119=8
label define visit_cross 2 "Month 1" 3 "Month 3" 4 "Month 6" 5 "Year 1" 6 "Year 2" 7 "Year 3" 8 "Year 8-9", replace
label values visit_cross visit_cross
tab visit_cross
tab follow_up_weeks_integer visit_cross if sampletype==1

replace visit = visit_cross if sampletype==1
drop visit_cross

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

//////////////////////////////////////////////////////////////////////////
*Merge in EE data
use "AIM 1 dataset_merged_eedata_December102025.dta", clear //Longitudinal dataset with EE matched-in
recode EE_score 0/6=0 7/10=1, gen(EE_bin)
tab EE_bin
summ concentration_NPO concentration_MPO concentration_AAT, d
keep pid visit EE_bin EE_score concentration_NPO concentration_MPO concentration_AAT
keep if visit==1
save "EE_lookup", replace
use "Both_datasets_long", clear
merge m:1 pid using "EE_lookup"
tab sampletype _merge, m
tab pid if _merge==2 //The EE dataset includes POL 421 who was later excluded
drop if _merge==2 
drop _merge

save "Both_datasets_long", replace
