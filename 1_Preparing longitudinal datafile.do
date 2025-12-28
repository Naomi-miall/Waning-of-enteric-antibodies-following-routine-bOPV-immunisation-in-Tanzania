/* Duration of polio antibody responses following bOPV vaccination in Tanzania.
  Alexandra Codeiro and Naomi Miall
  September 2025
 
*/

*SET UP
clear
set more off

*Filepath for Naomi:
cd "C:\Users\NaomiMiall\Documents\Polio work\Tanzania"

//////////////////////////////////////////////////////////////////
**# Anthro measures
*Anthro measures:
/* Indicate to the Stata compiler where the igrowup_restricted.ado file is stored*/
adopath + "Stata code\WHO anthro macro/"

/* Load the data file */
use poliovaccinestudy_Tanzania_aim1_baselinedata_092425, clear
describe

/* Generate the first three parameters reflib, datalib & datalab	*/

	//reflib - specify the file where the five STATA data sets containing the WHO Child Growth Standards are stored
	gen reflib="Stata code\WHO anthro macro/"
	lab var reflib "Directory of reference tables"

	//datalib - specify the file where the input STATA data set is stored.
	gen datalib="Stata code\Alex's code/"
	lab var datalib "Directory for datafiles"

	//datalab - specify the name that will prefix the output files (datalab_z_rc and datalab_prev_rc).
	gen datalab="Anthrooutputs_longitudinal"
	lab var datalab "Working file"

//For all variables ensure any missing are coded to ==.
/*	check the variable for "sex"	1 = male, 2=female */
gen sex = .
replace sex = 1 if femalesex == 0 // Can I double check coding of childgender? femalesex variable makes it look like 0=male?
replace sex = 2 if femalesex == 1  
label define sexlbl 1 "Male" 2 "Female"
label values sex sexlbl
tab sex femalesex, m

/*	check the variable for "age" */
**need to convert age in days to compare norms 
	gen agedays = interview_date - infant_dob
	label var agedays "Age in days at time of interview"
	summ agedays, d
	*checking:
		gen ageweek_integer = int(ageweek)
		tab agedays ageweek_integer

	/*	define your ageunit	*/
	gen ageunit="days"
	lab var ageunit "=days or =months"

/*	check the variable for body "weight" which must be in kilograms*/
/* 	NOTE: if not available, please create as [gen weight=.]*/
	summ weight_kgs, d

/* 	check the variable for "height" which must be in centimeters	*/ 
/* 	NOTE: if not available, please create as [gen height=.]*/
	summ height_cm, d

/*	check the variable for "measure"*/
//specify the name of a character variable indicating whether recumbent length or standing height was measured. The values of this variable must be "l" or "L" for recumbent length, and "h" or "H" for standing height. Users must code any missing values as " "
/* 	NOTE: if not available, please create as [gen str1 measure=" "]*/
	gen measure = ""
	replace measure ="l" if height_measurement_method==1
	replace measure = "h" if height_measurement_method == 2
	tab measure height_measurement_method, m
	label var measure "Length (l) or height (h) measured?"

//The macro will use height, measure and age to create _clenhei - the adjusted length/height with 0.7 where relevant.

/* 	check the variable for "oedema"*/
//values of this variable must be "n" or "N" for non-oedema, and "y" or "Y" for oedema. Macro assumes those with missing values are non-oedema. For oedema cases, weight-related z-scores (_zwei, _zwfl and _zbmi) are not calculated.
/* 	NOTE: if not available, please create as [gen str1 oedema="n"]*/
	tab legedema, m
	tab legedema, nolab //slight oedema == 1; none == 0
	gen oedema = ""
	replace oedema = "n" if legedema==0
	replace oedema = "y" if legedema==1
	tab legedema oedema, m

/*	check the variable for "sw" for the sampling weight*/
/* 	NOTE: if not available, please create as [gen sw=1] */
	gen sw=1
	
/* 	Fill in the macro parameters to run the command */
igrowup_restricted reflib datalib datalab sex agedays ageunit weight_kgs height_cm measure oedema sw

use "Stata code\Alex's code\Anthrooutputs_longitudinal_z_rc.dta", clear

*check height adjustment worked
*recommended that lenght measured for children < 731 days, if height measured for children < 731 days then add 0.7 cm. 
	summ agedays //all participants are < 731 days
	tab measure // had height measured so assume adjustment applied for 11
	list pid agedays measure height_cm _clenhei if height_cm!=_clenhei // looks good

*checking z scores worked
list pid agedays weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl in 1/15 //looks good

/* exploring extremes
_fwei=1 if (_zwei < -6 | _zwei >5) - extreme weight for age
_flen=1 if (_zlen < -6 | _zlen >6) - extreme height for age
_fbmi=1 if (_zbmi < -5 | _zbmi >5) - extreme bmi (for age)
_fwfl=1 if (_zwfl < -5 | _zwfl >5) - extreme weight for height
*/

list pid agedays weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if legedema==1

	summ _fwei _flen _fbmi _fwfl
	
	summ height_cm, d //median amongst participants is 61.5cm
	summ weight_kgs, d //median amongst participants is 6.6kgs
	
list pid agedays weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if _flen==1
list pid agedays weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if _fbmi==1
list pid agedays weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if _fwfl==1

/*Generating binary variables
WHO Indicator Definitions: 
	Stunting Height-for-age < -2SD
	Severe wasting Weight for Height < -3 SD
	Wasting Weight for Height < -2 SD
	Overweight Weight for Height > +2 SD
	Underweight Weight-for-age < -2SD
 SOURCE:https://cdn.who.int/media/docs/default-source/child-growth/child-growth-standards/software/anthro-survey-analyser-quickguide.pdf?sfvrsn=dc7ddc6f_6
*/

gen stunting = (_zlen<-2) if !missing(_zlen)
	label var stunting "Stunting defined as height-for-age below -2SD"
	tab stunting
	list pid sex agedays _clenhei _zlen if stunting==1

gen wasting = (_zwfl<-2) if !missing(_zwfl)
	label var wasting "Wasting defined as weight-for-height below -2SD"
	tab wasting
	list pid sex agedays height_cm weight_kgs _zwfl if wasting==1

gen severe_wasting = (_zwfl<-3) if !missing(_zwfl)
	label var severe_wasting "Severe wasting defined as weight-for-height below -3SD"
	tab severe_wasting
	list pid sex agedays height_cm weight_kgs _zwfl wasting if severe_wasting==1
	
gen underweight = (_zwei<-2) if !missing(_zwei)
	label var underweight "Underweight defined as weight-for-age below -2SD"
	tab underweight
	list pid sex agedays weight_kgs _zwei if underweight==1	

gen overweight = (_zwei>2) if !missing(_zwei)
	label var overweight "Overweight defined as weight-for-age above 2SD"
	tab overweight
	list pid sex agedays weight_kgs _zwei if overweight==1	
	
	drop sw measure oedema ageunit ageweek_integer agedays sex datalab datalib reflib _agedays _clenhei _cbmi _fwfl _flen _fwei _fbmi
	drop _zwei _zlen _zbmi
	
save "Stata code\Alex's code\Anthrooutputs_longitudinal_z_rc.dta", replace
//////////////////////////////////////////////////////////////////
**# Data cleaning

	*Creating co-infection exposure groupings and dropping repeat variables
use "Stata code\Alex's code\Anthrooutputs_longitudinal_z_rc.dta", clear

*HIV
tab hiv_pregnancy hiv_status_mother,m //hiv_pregnancy is repeat so keep only 1
drop hiv_pregnancy
rename hiv_status_mother hiv_exposure

*Amoeba
tab kato_katz_results parasite_identified_from_kato, m //repeats so keep only 1 and rename to coinfect_parasite for easy searching
rename kato_katz_results coinfect_parasite
label variable coinfect_parasite "Amoeba infection"
tab coinfect_parasite, m
drop parasite_identified_from_kato

*Schisto
tab serum_caa_results,m  // 2 missing, 14 positive
summ caa_concentration, d //large jump in values for the top 5%.
bysort serum_caa_results: summ caa_concentration //threshold for positive seems to be 30
//rename binary variable to coinfect_caa for easy seraching
rename serum_caa_results coinfect_caa
label variable coinfect_caa "Schisto infection (caa concentration > 30 in serum)"

*Group 1 - Enteroinvasive pathogens associated with dysentery (Enteroinvasive E. coli, Shigella-toxin producing E. coli, Campylobacter, Salmonella)
tab1 Enteroinvasive_ecoli_EIEC Shigella_toxin_STEC Campylobacter Salmonella //0; 2; 6; and 6 positive respectively. No missing data.
gen coinfect_enteroinvasive = 0
replace coinfect_enteroinvasive = 1 if Enteroinvasive_ecoli_EIEC==1 | Shigella_toxin_STEC==1 | Campylobacter==1 | Salmonella==1
label variable coinfect_enteroinvasive "Enteroinvasive pathogens associated with dysentery (Enteroinvasive E. coli, Shigella-toxin producing E. coli, Campylobacter, Salmonella)"
label values coinfect_enteroinvasive serum_caa
tab coinfect_enteroinvasive // 12 positive

*Group 2 - Other bacterial pathogens (C. difficile, P. shigelloides, Enterotoxigenic E. coli, Enteropathogenic E. coli, Enteroaggregative E. coli)
tab1 Clostridium_difficile Plesiomonas_shigelloides Enterotoxigenic_ecoli_ETEC Enteropathogenic_ecoli_EPEC Enteroaggregative_ecoli_EAEC //3; 3; 11; 20 and 69 positive respectively. No missing data.
gen coinfect_otherbacteria = 0
replace coinfect_otherbacteria = 1 if Clostridium_difficile==1 | Plesiomonas_shigelloides==1 | Enterotoxigenic_ecoli_ETEC==1 | Enteropathogenic_ecoli_EPEC==1 | Enteroaggregative_ecoli_EAEC==1 
label variable coinfect_otherbacteria "Other bacterial pathogens (C. difficile, P. shigelloides, Enterotoxigenic E. coli, Enteropathogenic E. coli, Enteroaggregative E. coli)"
label values coinfect_otherbacteria serum_caa
tab coinfect_otherbacteria // 73 positive

*Group 3 - Rotavirus, Adenovirus, Sapovirus, Astrovirus, or Norovirus
tab1 Rotavirus_A Adenovirus Astrovirus Norovirus Sapovirus // 27; 2; 3; 13 and 2 positive respectively. No missing data.
gen coinfect_virus = 0
replace coinfect_virus = 1 if Rotavirus_A==1 | Adenovirus==1 | Astrovirus==1 | Norovirus==1 | Sapovirus==1
label variable coinfect_virus "Rotavirus, Adenovirus, Sapovirus, Astrovirus, or Norovirus"
label values coinfect_virus serum_caa
tab coinfect_virus // 41 positive - minimal overlapp in these viruses

*Giardia lamblia
rename Giardia_lamblia coinfect_giardia
label variable coinfect_giardia "Giardia lamblia"
label values coinfect_giardia serum_caa
tab coinfect_giardia // 5 positive

*Other
tab1 Cryptosporidium Cyclospora_cayetanensis Entamoeba_histolytica Vibrio Yersinia_enterocolitica //No participants have any of the above


*Any infection and number of infections (including parasites, schisto, biofire)
*NB - existing variable None_detected refers to biofire, not to schisto
tab1 coinfect_caa coinfect_enteroinvasive coinfect_giardia coinfect_otherbacteria coinfect_parasite coinfect_virus if None_detected==1
egen no_infections_ungrouped =  anycount(coinfect_parasite coinfect_caa coinfect_giardia coinfect_enteroinvasive coinfect_otherbacteria coinfect_virus), values(1)
tab no_infections_ungrouped
label variable no_infections_ungrouped "Number of infections ungrouped, not including HIV exposure"
gen any_infection = no_infections_ungrouped
replace any_infection = 1 if no_infections_ungrouped>0

drop Enteroinvasive_ecoli_EIEC Shigella_toxin_STEC Campylobacter Salmonella Clostridium_difficile Plesiomonas_shigelloides Enterotoxigenic_ecoli_ETEC Enteropathogenic_ecoli_EPEC Enteroaggregative_ecoli_EAEC Adenovirus Astrovirus Norovirus Sapovirus Cryptosporidium Cyclospora_cayetanensis Entamoeba_histolytica Vibrio Yersinia_enterocolitica

*Cohort
gen sampletype=2
label define sampletype 1 "cross-sectional" 2 "longitudinal"
label values sampletype sampletype

save "Baseline_data_preppedformerge", replace
use "Baseline_data_preppedformerge", clear
/////////////////////////////////////////////////////////
**# Merging neuts data
*Study samples data:
*1. Prepping dates file in long
*N.b. Make sure no space after POL 540 in the excel file
import excel "Stata code\Polio vaccine study follow up dates.xlsx", sheet("Sheet1") cellrange(A1:D104) firstrow case(lower) clear

rename monthfollowupdateofvisit date2
rename c date3
rename d date4

reshape long date, i(pid) j(visit)

save "follow_up_dates", replace

*prepare neuts data file for merge
use "stool_serumneutralization_09242025.dta", clear

*Make variables for number of observations, visit number, and binary variables for attendance at each follow-up
tab study_visitmonth, m //103 at enrollment, 84 at 1 month, 79 at 3, 77 at 6
bysort pid: egen no_obs = count(pid)
gen visit = 1
replace visit = 2 if study_visitmonth== "1 month"
replace visit = 3 if study_visitmonth== "3 months"
replace visit = 4 if study_visitmonth== "6 months"
//Make a binary variable for which month each obs is in
gen one_month_p = 0
replace one_month_p = 1 if study_visitmonth== "1 month"
bysort pid: egen one_month = max(one_month_p)
label variable one_month "observed (non-missing) at 1 month"
gen three_months_p = 0
replace three_months_p = 1 if study_visitmonth== "3 months"
bysort pid: egen three_months = max(three_months_p)
label variable three_months "observed (non-missing) at 3 months"
gen six_months_p = 0
replace six_months_p  = 1 if study_visitmonth== "6 months"
bysort pid: egen six_months = max(six_months_p)
label variable six_months "observed (non-missing) at 6 months"
drop study_visitmonth one_month_p three_months_p six_months_p

*merge in long format
merge 1:1 pid visit using "follow_up_dates"
sort pid visit 
drop if no_obs==.
list pid visit date one_month three_months six_months no_obs in 1/50, sepby(pid)

list pid if _merge==1 & visit>1 //none
drop _merge

tab no_obs if visit==1 //66 fully complete
tab1 one_month three_months six_months if visit==1 //103 in visit 1, 84 in visit 2, 79 in 3, 77 in 4.

merge m:1 pid using "Baseline_data_preppedformerge", nogen

replace date = interview_date if date==. & visit==1
drop interview_date
list in 1

*Data cleaning of mucosal immunity data:
**convert titers to log(2) scale. Higher titer indicates higher amount of antibodies. The oral poliovirus vaccine generates antibodies to PV1 and PV3 (types 1 and 3 poliovirus)
**Using the equation log2(x) = ln(x)/ln(2)
generate logPV1_stool = ln(pv1_stool)/ln(2)
generate logPV2_stool = ln(pv2_stool)/ln(2)
generate logPV3_stool = ln(pv3_stool)/ln(2)

generate logPV1_serum = ln(pv1neuttiterserum)/ln(2)
generate logPV2_serum = ln(pv2neuttiterserum)/ln(2)
generate logPV3_serum = ln(pv3neuttiterserum)/ln(2)

*Generating binary neuts data:
foreach x in 1 2 3 {
	gen pv`x'_stool_positive = 0
	replace pv`x'_stool_positive = 1 if pv`x'_stool>=8
	replace pv`x'_stool_positive = . if pv`x'_stool==.
	bysort pv`x'_stool_positive: summ pv`x'_stool 
	tab pv`x'_stool_positive if pv`x'_stool==., m
	gen pv`x'_stool_detectable = 0
	replace pv`x'_stool_detectable = 1 if pv`x'_stool>=4
	replace pv`x'_stool_detectable = . if pv`x'_stool==.
	bysort pv`x'_stool_detectable: summ pv`x'_stool
}
tab1 pv1_stool_positive pv2_stool_positive pv3_stool_positive pv1_stool_detectable pv2_stool_detectable pv3_stool_detectable
tab pv1_stool_positive pv3_stool_positive, col

*Make a strata variable - high or low at baseline for neuts
foreach x in 1 2 3 {
	gen high_neut_baseline_pv`x'_prep = pv`x'_stool_positive
	replace high_neut_baseline_pv`x'_prep = . if visit>1
	bysort pid: egen high_neut_baseline_pv`x' = max(high_neut_baseline_pv`x'_prep)
	drop high_neut_baseline_pv`x'_prep
	sort pid visit
	list pid visit pv`x'_stool_positive high_neut_baseline_pv`x' in 1/20
}

//A few missing:
list pid visit if pv1_stool==. //4 obs missing.

//Tag if missing any stool neut data
tab logPV1_stool if logPV3_stool == ., m
tab logPV3_stool if logPV1_stool == ., m
tab logPV2_stool if logPV3_stool == ., m //same particpants are missing on all
bysort pid: egen flag_missing_stool_neuts = max(logPV1_stool), missing
recode flag_missing_stool_neuts 0/10=0 .=1
sort pid visit
list visit flag_missing_stool_neuts logPV1_stool logPV3_stool logPV3_stool if flag_missing_stool_neuts==1, sepby(pid) //just missing on 1 obs, not all of them.

//Tag if missing any serum neut data
 misstable summarize logPV1_serum logPV2_serum logPV3_serum
tab logPV1_serum if logPV3_serum == ., m
tab logPV3_serum if logPV1_serum == ., m
tab logPV2_serum if logPV3_serum == ., m
 //again same particpants are missing on all
tab visit if logPV3_serum == ., m // spread across all visits
tab logPV3_serum if logPV3_stool==., m //3 of the obs missing stool neut are also missing serum, 1 has valid serum.
 
bysort pid: egen flag_missing_serum_neuts = max(logPV1_serum), missing
recode flag_missing_serum_neuts 0/50=0 .=1
sort pid visit
list visit flag_missing_serum_neuts logPV1_serum logPV3_serum logPV3_serum if flag_missing_serum_neuts==1, sepby(pid)

//Binary - response to PV1 or PV3 at any time-point inlcuding baseline
gen pv1_stool_responses_prep = 0
replace pv1_stool_responses_prep = . if logPV1_stool==.
replace pv1_stool_responses_prep = 1 if pv1_stool_positive==1
bysort pid: egen pv1_stool_responses = max(pv1_stool_responses_prep)
drop pv1_stool_responses_prep 
tab pv1_stool_responses  if visit==1

gen pv1_stool_responses_exclbl_prep = 0
replace pv1_stool_responses_exclbl_prep = . if logPV1_stool==. & visit>1
replace pv1_stool_responses_exclbl_prep = 1 if pv1_stool_positive==1 & visit>1
bysort pid: egen pv1_stool_responses_exclbl = max(pv1_stool_responses_exclbl_prep)
drop pv1_stool_responses_exclbl_prep 
tab pv1_stool_responses_exclbl  if visit==1

gen pv3_stool_responses_prep = 0
replace pv3_stool_responses_prep = . if logPV3_stool==.
replace pv3_stool_responses_prep = 1 if pv3_stool_positive==1
bysort pid: egen pv3_stool_responses = max(pv3_stool_responses_prep)
drop pv3_stool_responses_prep 
tab pv3_stool_responses  if visit==1

gen pv3_stool_responses_exclbl_prep = 0
replace pv3_stool_responses_exclbl_prep = . if logPV3_stool==. & visit>1
replace pv3_stool_responses_exclbl_prep = 1 if pv3_stool_positive==1 & visit>1
bysort pid: egen pv3_stool_responses_exclbl = max(pv3_stool_responses_exclbl_prep)
drop pv3_stool_responses_exclbl_prep 
tab pv3_stool_responses_exclbl  if visit==1


save "merged_dataset_long", replace
use "merged_dataset_long", clear


/////////////////////////////////////////////////////////
*Prepping variables for Table 1
*Mother age at birth
gen age_mother_weeks = (age_mother*52)+26
gen age_mother_at_birth = int((age_mother_weeks - ageweek)/52)
list age_mother_at_birth ageweek age_mother in 1/40, sepby(pid)

*Cat birthweight
summarize birthweight, detail
describe birthweight
destring birthweight, replace
gen bw_bucket = .
replace bw_bucket = 1 if birthweight < 1.5        
replace bw_bucket = 2 if birthweight >= 1.5 & birthweight < 2.5 
replace bw_bucket = 3 if birthweight >= 2.5
label define bw_cat 1 "Very low (<1.5)" 2 "<2.5" 3 "Normal (≥2.5)"
label values bw_bucket bw_cat
tab bw_bucket
tabulate vaccines_not_received 

*bOPV
gen age_bopv0_days = opv_date - infant_dob
gen delay_bopv0 = .
replace delay_bopv0 = 1 if age_bopv0_days > 7
replace delay_bopv0 = 0 if age_bopv0_days <= 7 & !missing(age_bopv0_days)
label define delay_lbl 0 "On time (≤7d)" 1 "Delayed (>7d)"
label values delay_bopv0 delay_lbl
tabulate delay_bopv0, m

**calculating age in days at the first bOPV vaccination 
gen age_opv0_days = opv_date - infant_dob
tabulate age_opv0_days
summarize age_opv0_days, detail 

**weeks from birth till bOPV 4 at 14 weeks
**calculating age in weeks at the time of week 14 vaccination
gen age_opv14_days = bopv14w_date - infant_dob
gen age_opv14_weeks = round(age_opv14_days / 7, 0.1)
tabulate age_opv14_weeks
summarize age_opv14_weeks, detail

**calculating age in weeks at the time of IPV vaccine
gen ipv_days = ipv_date - infant_dob
gen ipv_weeks = round(ipv_days / 7, 0.1)
tabulate ipv_weeks
summarize ipv_weeks, detail 

drop opv_date bopv6w_date bopv10w_date ipv_date ipv_days age_opv0_days age_opv6_days age_opv10_days age_opv14_days

**ROTAVIRUS data
**Looking at the days/months since rotavirus vaccine. 
gen days_since_rota6 = date - rota_vaccine_w6
gen days_since_rota10 = date - rota_vaccine_w10
gen weeks_since_rota10 = int(days_since_rota10/7)
summarize days_since_rota6 days_since_rota10, detail

*Make variable for follow-up length in days
gen follow_up_days = .
replace follow_up_days = date - bopv14w_date
list bopv14w_date date follow_up_days age_opv14_weeks in 1/10
histogram follow_up_days, bin(200)

*Flagging age outliers
tab age_opv14_weeks if visit==1
hist ageweek if visit==1, bin(100)
gen flag_ageoutlier = 0
replace flag_ageoutlier = 1 if age_opv14_weeks > 40 

*variable for average 'day' for each visit (needed for summary stats on plots)
bysort visit: summ follow_up_days
gen summ_visit_day = .
replace summ_visit_day = 0 if visit==1
replace summ_visit_day = 30 if visit==2
replace summ_visit_day = 90 if visit==3
replace summ_visit_day = 180 if visit==4

//Binary variables to simplify regression analyses
gen ageweek_bin = 0
replace ageweek_bin = 1 if ageweek>=17.14286 //upper quartile
label define ageweek_binl 0 "<17.14 weeks" 1 ">=17.14 weeks, Upper quartile", replace
label values ageweek_bin ageweek_binl
tab ageweek_bin, m

gen notbreastfeeding_bin = infant_feeding_type
recode notbreastfeeding_bin 1=0 2=1 4=1
label define notbreastfeeding_binl 0 "Exclusive breastfeeding" 1 "Not exclusive"
label values notbreastfeeding_bin notbreastfeeding_binl
tab notbreastfeeding_bin

gen bf_4months=notbreastfeeding_bin
recode bf_4months 0=1 1=0

gen watersource_bin = 0
replace watersource_bin = 1 if minutes_towatersource>=15 //upper quartile
label define watersource_binl 0 "<15 minutes" 1 ">=15 minutes, Upper quartile"
label values watersource_bin watersource_binl
tab watersource_bin

gen summ_visit_month = .
replace summ_visit_month = 0 if visit==1
replace summ_visit_month = 1 if visit==2
replace summ_visit_month = 3 if visit==3
replace summ_visit_month = 6 if visit==4

gen follow_up_months = follow_up_days/30 

save "merged_dataset_long", replace
use "merged_dataset_long", clear

//////////////////////////////////////////////
**# Merging in Ig data
*Preparing Stool IgA datafile
	*Tidy up pid values and add sampletype and visit labels
use "stool_IgA_withlabids_studyvisitmonths.dta", clear
replace pid = "POL 482" if pid=="POL. 482"
replace pid = "POL 282" if pid=="POL 282 "
drop if pid=="POL 421" //not included in the baseline results / had been dropped from the study
tab study_visit_month
gen sampletype = 2
replace sampletype = 1 if study_visit_month=="aim2"
label define sampletype 1 "cross-sectional" 2 "longitudinal"
label values sampletype sampletype
tab sampletype
keep if sampletype==2
gen visit = 1
replace visit = 2 if study_visit_month=="1-month"
replace visit = 3 if study_visit_month=="3-months"
replace visit = 4 if study_visit_month=="6-month" 
drop study_visit_month
duplicates tag pid visit, gen(duplicated_flag)
sort pid visit
list if duplicated>0, sepby(pid)
drop if duplicated==3

	*For samples where IgA was measured twice for same individual and timepoint, take the mean and then transform all to log10 scale.
foreach x in 1 2 3{
	destring IPV`x', replace
	bysort pid visit: egen PV`x'_IgA_stool = mean(IPV`x')
	list in 1/10, sepby(pid)
	list if duplicated>0, sepby(pid)
	drop IPV`x'
	gen logPV`x'_IgA_stool = ln(PV`x'_IgA_stool)/ln(10)
}
list if duplicated>0, sepby(pid)
duplicates drop pid visit, force
list in 1/10, sepby(pid)
save "longitudinal_stool_IgA", replace

*preparing Serum IgG file
use "serum_IgG_1000f_with_labids_studyvisits.dta", clear //Using 1000f dilution rather than 10,000f due to better sensivity (less saturation at the top)
replace pid = "POL 482" if pid=="POL. 482"
tab study_visit_month
gen sampletype = 2
replace sampletype = 1 if study_visit_month=="Aim2"
label define sampletype 1 "cross-sectional" 2 "longitudinal"
label values sampletype sampletype
tab sampletype
keep if sampletype==2
gen visit = 1
replace visit = 2 if study_visit_month=="1-month"
replace visit = 3 if study_visit_month=="3-month "
replace visit = 4 if study_visit_month=="6-months" 
tab study_visit_month visit, m
drop study_visit_month
duplicates tag pid visit, gen(duplicated_flag)
sort pid visit
list if duplicated>0, sepby(pid)

foreach x in 1 2 3{
	destring IPV`x', replace
	bysort pid visit: egen PV`x'_IgG_serum  = mean(IPV`x')
	list in 1/10, sepby(pid)
	list if duplicated>0, sepby(pid)
	drop IPV`x'
	gen logPV`x'_IgG_serum = ln(PV`x'_IgG_serum )/ln(10)
}
list if duplicated>0, sepby(pid)
duplicates drop pid visit, force
list in 1/10, sepby(pid)
replace pid = "POL 482" if pid=="POL. 482"
save "longitudinal_serum_IgG", replace

*preparing Serum IgA file
use "serum_IgA_with_labids_studyvisits.dta", clear
replace pid = "POL 482" if pid=="POL. 482"
tab study_visit_month
gen sampletype = 2
replace sampletype = 1 if study_visit_month=="Aim2"
label define sampletype 1 "cross-sectional" 2 "longitudinal"
label values sampletype sampletype
tab sampletype
keep if sampletype==2
gen visit = 1
replace visit = 2 if study_visit_month=="1-month"
replace visit = 3 if study_visit_month=="3-months"
replace visit = 4 if study_visit_month=="6-months" 
tab study_visit_month visit, m
drop study_visit_month
duplicates tag pid visit, gen(duplicated_flag)
sort pid visit
list if duplicated>0, sepby(pid)

foreach x in 1 2 3{
	replace IPV`x'="." if IPV`x'=="NaN"
	destring IPV`x', replace
	bysort pid visit: egen PV`x'_IgA_serum  = mean(IPV`x')
	list in 1/10, sepby(pid)
	list if duplicated>0, sepby(pid)
	drop IPV`x'
	gen logPV`x'_IgA_serum = ln(PV`x'_IgA_serum)/ln(10)
}
list if duplicated>0, sepby(pid)
duplicates drop pid visit, force
list in 1/10, sepby(pid)

save "longitudinal_serum_IgA", replace

*Merge all IgA and IgG data
use "longitudinal_serum_IgA", clear
merge 1:1 pid visit using "longitudinal_serum_IgG.dta"
list pid sampletype visit duplicate logPV1* if _merge!=3
drop _merge

merge 1:1 pid visit using "longitudinal_stool_IgA.dta"
sort pid visit
list pid sampletype visit duplicate logPV1* if _merge!=3, sepby(pid)
drop _merge
save "longitudinal_Ig_data", replace
use "longitudinal_Ig_data", clear
keep if visit==1
duplicates report pid //103 pids - correct.

use "merged_dataset_long", clear
merge 1:1 pid visit using "longitudinal_Ig_data"
list pid sampletype visit duplicate logPV1* if _merge!=3, sepby(pid)
drop _merge

save "merged_dataset_long", replace

